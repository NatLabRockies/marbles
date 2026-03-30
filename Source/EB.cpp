#include "EB.H"
#include "Geometry.H"
#include <cstring> // for std::memcpy
#include <sstream> // for CSV parsing

namespace lbm {
void initialize_eb(const amrex::Geometry& geom, const int max_level)
{
    BL_PROFILE("LBM::initialize_eb()");

    amrex::ParmParse pp("eb2");

    std::string geom_type("all_regular");
    pp.query("geom_type", geom_type);

    int max_coarsening_level = max_level;
    amrex::ParmParse ppamr("amr");
    amrex::Vector<int> ref_ratio(max_level, 2);
    ppamr.queryarr("ref_ratio", ref_ratio, 0, max_level);
    for (int lev = 0; lev < max_level; ++lev) {
        max_coarsening_level +=
            (ref_ratio[lev] == 2
                 ? 1
                 : 2); // Since EB always coarsening by factor of 2
    }

    // Custom types defined here - all_regular, plane, sphere, etc, will get
    // picked up by default (see AMReX_EB2.cpp around L100 )
    amrex::Vector<std::string> amrex_defaults(
        {"all_regular", "box", "cylinder", "plane", "sphere", "torus", "parser",
         "stl"});
    if (std::find(amrex_defaults.begin(), amrex_defaults.end(), geom_type) ==
        amrex_defaults.end()) {
        std::unique_ptr<lbm::Geometry> geometry(
            lbm::Geometry::create(geom_type));
        geometry->build(geom, max_coarsening_level);
    } else {
        // For all AMReX default types (including voxel_cracks), use standard
        // build voxel_cracks will override the m_is_fluid in
        // initialize_from_stl
        amrex::EB2::Build(geom, max_level, max_level);
    }
}

void initialize_from_stl(
    const amrex::Geometry& geom, amrex::iMultiFab& is_fluid)
{
    BL_PROFILE("LBM::initialize_from_stl()");

    amrex::ParmParse pp("eb2");
    std::string geom_type("all_regular");
    pp.query("geom_type", geom_type);
    std::string name;
    pp.query("stl_file", name);

    // use native AMReX EB STL utility
    if ((!name.empty()) && (geom_type == "stl")) {
        return;
    }

    if ((!name.empty()) && (geom_type == "all_regular")) {
        amrex::Real scale = 1.0;
        int reverse_normal = 0;
        amrex::Array<amrex::Real, 3> center = {0.0, 0.0, 0.0};
        pp.query("stl_scale", scale);
        pp.query("stl_reverse_normal", reverse_normal);
        pp.query("stl_center", center);

        amrex::STLtools stlobj;
        stlobj.read_stl_file(name, scale, center, reverse_normal);

        amrex::MultiFab marker(
            is_fluid.boxArray(), is_fluid.DistributionMap(), 1,
            is_fluid.nGrow());

        const amrex::Real outside_value = 1.0;
        const amrex::Real inside_value = 0.0;
        marker.setVal(1.0);
        stlobj.fill(
            marker, marker.nGrowVect(), geom, outside_value, inside_value);
        amrex::Gpu::synchronize();

        auto const& marker_arrs = marker.const_arrays();
        auto const& is_fluid_arrs = is_fluid.arrays();
        amrex::ParallelFor(
            is_fluid, is_fluid.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                is_fluid_arrs[nbx](i, j, k, 0) =
                    static_cast<int>(marker_arrs[nbx](i, j, k, 0));
            });
        amrex::Gpu::synchronize();
    }

    // Check for voxel crack generation flag
    int use_voxel_cracks = 0;
    pp.query("use_voxel_cracks", use_voxel_cracks);
    if (use_voxel_cracks != 0) {
        amrex::Print() << "Using voxel crack generation" << std::endl;
        generate_voxel_cracks(geom, is_fluid);
        return;
    }

    if ((!name.empty()) && (geom_type != "all_regular")) {
        amrex::Abort(
            "LBM::initialize_from_stl() geom_type should be all_regular to "
            "avoid issues");
    }
}

std::vector<uint16_t>
read_crack_file(const std::string& filename, int nx, int ny, int nz)
{
    BL_PROFILE("LBM::read_crack_file()");

    // Auto-detect file format
    bool is_csv =
        (filename.size() >= 4 &&
         filename.substr(filename.size() - 4) == ".csv");

    // Use AMReX's cross-platform file reading utilities
    amrex::Vector<char> file_char_ptr;

    // Use AMReX's parallel-safe file reading
    try {
        amrex::ParallelDescriptor::ReadAndBcastFile(filename, file_char_ptr);
    } catch (const std::exception& e) {
        amrex::Abort(
            "Error reading crack file: " + filename + " - " +
            std::string(e.what()));
    }

    std::vector<uint16_t> crack_data(
        static_cast<size_t>(nx) * static_cast<size_t>(ny) *
        static_cast<size_t>(nz));

    if (is_csv) {
        // Parse CSV format
        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::fill(
                crack_data.begin(), crack_data.end(), 1); // Initialize as solid

            std::string file_content(
                file_char_ptr.data(), file_char_ptr.size());
            std::istringstream iss(file_content);
            std::string line;

            // Skip header line
            if (!std::getline(iss, line)) {
                amrex::Abort("CSV file is empty or corrupted: " + filename);
            }

            size_t line_count = 0;
            while (std::getline(iss, line)) {
                // Skip empty lines
                if (line.empty()) {
                    continue;
                }

                std::istringstream line_stream(line);
                std::string token;

                // Parse tokens: X, Y, Z (discard), then tag
                int tag;
                try {
                    // discard X
                    if (!std::getline(line_stream, token, ',')) {
                        continue;
                    }
                    // discard Y
                    if (!std::getline(line_stream, token, ',')) {
                        continue;
                    }
                    // discard Z
                    if (!std::getline(line_stream, token, ',')) {
                        continue;
                    }
                    // read tag
                    if (std::getline(line_stream, token, ',')) {
                        tag = std::stoi(token);
                    } else {
                        continue;
                    }
                } catch (const std::exception& e) {
                    // Skip malformed lines
                    continue;
                }

                // CSV and binary are written in identical sequence by
                // mainCrackGenerator.C Both loop: k(Z) -> j(Y) -> i(X), so just
                // read sequentially The X,Y,Z coordinates are metadata - what
                // matters is the order

                if (line_count < static_cast<size_t>(nx) *
                                     static_cast<size_t>(ny) *
                                     static_cast<size_t>(nz)) {
                    crack_data[line_count] = static_cast<uint16_t>(tag);
                }
                line_count++;
            }

            amrex::Print() << "Successfully read CSV crack file: " << filename
                           << " (" << line_count << " data points)"
                           << std::endl;
        }
        // Broadcast the parsed data to all processors
        amrex::ParallelDescriptor::Bcast(
            crack_data.data(), crack_data.size() * sizeof(uint16_t), 0);

    } else {
        // Parse binary format
        size_t expected_size = static_cast<size_t>(nx) *
                               static_cast<size_t>(ny) *
                               static_cast<size_t>(nz) * sizeof(uint16_t);

        // ReadAndBcastFile may add extra bytes, so ensure we only use what we
        // need
        if (static_cast<size_t>(file_char_ptr.size()) < expected_size) {
            amrex::Abort(
                "Binary file too small after reading. Expected: " +
                std::to_string(expected_size) + " bytes, got: " +
                std::to_string(file_char_ptr.size()) + " bytes");
        }

        // Convert char data to uint16_t array
        std::memcpy(crack_data.data(), file_char_ptr.data(), expected_size);

        if (amrex::ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "Successfully read binary crack file: "
                           << filename << " (" << expected_size << " bytes)"
                           << std::endl;
        }
    }

    return crack_data;
}

void generate_voxel_cracks(
    const amrex::Geometry& geom, amrex::iMultiFab& is_fluid)
{
    BL_PROFILE("LBM::generate_voxel_cracks()");

    amrex::ParmParse pp("voxel_cracks");

    // Get grid dimensions from domain
    const amrex::Box& domain = geom.Domain();
    const int nx = domain.length(0);
    const int ny = domain.length(1);
    const int nz = domain.length(2);

    amrex::Print() << "Loading voxel cracks for domain: " << nx << " x " << ny
                   << " x " << nz << std::endl;

    // Get binary crack file path
    std::string crack_file;
    if (pp.query("crack_file", crack_file) == 0) {
        // Default filename pattern matching mainCrackGenerator.C output
        crack_file = "microstructure_nX" + std::to_string(nx) + "_nY" +
                     std::to_string(ny) + "_nZ" + std::to_string(nz) + ".bin";
    }

    // Read crack pattern from file (auto-detect format)
    std::vector<uint16_t> crack_data = read_crack_file(crack_file, nx, ny, nz);

    // Initialize all cells as SOLID first
    is_fluid.setVal(1);

    // Copy crack data to GPU-accessible memory for CUDA builds
    amrex::Gpu::DeviceVector<uint16_t> d_crack_data(crack_data.size());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, crack_data.begin(), crack_data.end(),
        d_crack_data.begin());
    amrex::Gpu::synchronize();

    // Get raw pointer for device access
    auto const* crack_ptr = d_crack_data.data();

    // Copy crack data to MultiFab using GPU-compatible approach
    // Your file stores in k,j,i (z,y,x) order
    for (amrex::MFIter mfi(is_fluid); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox();
        amrex::Array4<int> const& is_fluid_arr = is_fluid.array(mfi);

        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Convert AMReX (i,j,k) to file index (z,y,x order)
                int file_index = k * (nx * ny) + j * nx + i;
                // Your binary file: 0 = fluid (tubes), 1 = solid
                // AMReX m_is_fluid: 0 = solid, 1 = fluid
                // So we need to invert the values
                is_fluid_arr(i, j, k, 0) = (crack_ptr[file_index] == 0) ? 1 : 0;
            });
    }
    amrex::Gpu::synchronize();

    amrex::Print() << "Voxel crack generation complete" << std::endl;
}

} // namespace lbm
