#include "EB.H"
#include "Geometry.H"
#include "Constants.H"
#include <cstring>   // for std::memcpy
#include <sstream>   // for CSV parsing
#include <fstream>   // for STL parsing in voxelize_stl_double_precision
#include <cstdint>   // for uint32_t in STL binary reader
#include <algorithm> // for std::min / std::max in triangle bbox reject
#include <array>     // for std::array in the double-precision STL reader

namespace lbm {

namespace {

// ---------------------------------------------------------------
// Double-precision STL reader + segment-triangle intersect helpers
// used by voxelize_stl_double_precision.  Kept in an anonymous
// namespace so they are TU-local and do not appear in EB.H.
// ---------------------------------------------------------------

using DTri = std::array<std::array<double, 3>, 3>; // 3 vertices, 3 coords each

// Read a binary STL file.  Binary STL format (little-endian):
//   80  bytes  header
//    4  bytes  uint32 triangle count
// per triangle:
//   12  bytes  normal   (3 x float32, discarded -- computed from vertices)
//   12  bytes  vertex 1 (3 x float32)
//   12  bytes  vertex 2 (3 x float32)
//   12  bytes  vertex 3 (3 x float32)
//    2  bytes  attribute byte count (discarded)
static bool
read_binary_stl_double(const std::string& fname, std::vector<DTri>& tris)
{
    std::ifstream is(fname, std::ios::binary);
    if (!is.good()) {
        return false;
    }
    char header[80];
    is.read(header, 80);
    uint32_t ntri = 0;
    is.read(reinterpret_cast<char*>(&ntri), 4);
    if (!is) {
        return false;
    }
    tris.resize(ntri);
    for (uint32_t t = 0; t < ntri; ++t) {
        float discard[3];
        is.read(reinterpret_cast<char*>(discard), 12); // normal
        for (int v = 0; v < 3; ++v) {
            float vf[3];
            is.read(reinterpret_cast<char*>(vf), 12);
            tris[t][v][0] = static_cast<double>(vf[0]);
            tris[t][v][1] = static_cast<double>(vf[1]);
            tris[t][v][2] = static_cast<double>(vf[2]);
        }
        char attr[2];
        is.read(attr, 2);
        if (!is) {
            return false;
        }
    }
    return true;
}

// Read an ASCII STL file.  Grammar:
//   solid <name>
//     facet normal <nx> <ny> <nz>
//       outer loop
//         vertex <x> <y> <z>
//         vertex <x> <y> <z>
//         vertex <x> <y> <z>
//       endloop
//     endfacet
//     ...
//   endsolid
static bool
read_ascii_stl_double(const std::string& fname, std::vector<DTri>& tris)
{
    std::ifstream is(fname);
    if (!is.good()) {
        return false;
    }
    std::string tok;
    is >> tok; // "solid"
    if (tok != "solid") {
        return false;
    }
    // Consume the rest of the header line (solid name).
    std::string rest;
    std::getline(is, rest);

    while (is >> tok) {
        if (tok == "endsolid") {
            break;
        }
        if (tok != "facet") {
            continue;
        }
        // facet normal nx ny nz
        is >> tok; // "normal"
        double nx, ny, nz;
        is >> nx >> ny >> nz;
        // outer loop
        is >> tok >> tok;
        DTri tri;
        for (int v = 0; v < 3; ++v) {
            is >> tok; // "vertex"
            is >> tri[v][0] >> tri[v][1] >> tri[v][2];
        }
        is >> tok; // "endloop"
        is >> tok; // "endfacet"
        tris.push_back(tri);
    }
    return true;
}

// Segment [a, b] vs triangle [t1, t2, t3] intersection in double.
//
// Algorithm: plane test for the segment endpoints, then a 3D
// barycentric-coordinate check for the intersection point.
// Tie-breaking convention matches amrex::STLtools: intersections
// exactly on triangle edges / vertices count once.
//
// Reference for the barycentric formulas:
//   Real-Time Collision Detection, C. Ericson (2005), sec. 5.1.5.
//
// Robustness: no epsilons are required at double precision for a
// randomly-chosen reference point far outside the mesh (the caller
// picks such a point in `voxelize_stl_double_precision`); any
// coincidental edge hit would be perturbed by O(1e-16) and land
// unambiguously on one side.
static bool segment_triangle_intersect_double(
    const double a[3], const double b[3], const DTri& tri)
{
    const double* t1 = tri[0].data();
    const double* t2 = tri[1].data();
    const double* t3 = tri[2].data();

    // Bounding-box quick reject (per-axis, both directions).
    for (int d = 0; d < 3; ++d) {
        const double tlo = std::min({t1[d], t2[d], t3[d]});
        const double thi = std::max({t1[d], t2[d], t3[d]});
        const double slo = std::min(a[d], b[d]);
        const double shi = std::max(a[d], b[d]);
        if (shi < tlo || slo > thi) {
            return false;
        }
    }

    // Edge vectors and triangle normal (unnormalized).
    const double e1[3] = {t2[0] - t1[0], t2[1] - t1[1], t2[2] - t1[2]};
    const double e2[3] = {t3[0] - t1[0], t3[1] - t1[1], t3[2] - t1[2]};
    const double n[3] = {
        e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2],
        e1[0] * e2[1] - e1[1] * e2[0]};

    // Segment direction.
    const double d[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};

    // Segment ~ parallel to triangle plane -> no unique intersection.
    const double dot_dn = d[0] * n[0] + d[1] * n[1] + d[2] * n[2];
    if (dot_dn == 0.0) {
        return false;
    }

    // Solve for parameter u in [0, 1] such that P = a + u*d lies on the
    // triangle's plane.  P - t1 . n = 0  =>  (a-t1 + u d) . n = 0.
    const double au_t1[3] = {t1[0] - a[0], t1[1] - a[1], t1[2] - a[2]};
    const double u_num = au_t1[0] * n[0] + au_t1[1] * n[1] + au_t1[2] * n[2];
    const double u = u_num / dot_dn;
    if (u < 0.0 || u > 1.0) {
        return false;
    }

    // Intersection point in triangle plane.
    const double P[3] = {a[0] + u * d[0], a[1] + u * d[1], a[2] + u * d[2]};

    // Barycentric coordinates via dot products (Ericson 2005).
    const double p1[3] = {P[0] - t1[0], P[1] - t1[1], P[2] - t1[2]};
    const double dot11 = e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2];
    const double dot12 = e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2];
    const double dot22 = e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2];
    const double dot1P = e1[0] * p1[0] + e1[1] * p1[1] + e1[2] * p1[2];
    const double dot2P = e2[0] * p1[0] + e2[1] * p1[1] + e2[2] * p1[2];
    const double denom = dot11 * dot22 - dot12 * dot12;
    if (denom == 0.0) {
        return false; // degenerate triangle
    }
    const double inv = 1.0 / denom;
    const double bary1 = (dot22 * dot1P - dot12 * dot2P) * inv;
    const double bary2 = (dot11 * dot2P - dot12 * dot1P) * inv;

    return (bary1 >= 0.0) && (bary2 >= 0.0) && (bary1 + bary2 <= 1.0);
}

} // anonymous namespace

// ---------------------------------------------------------------
// Public entry point (see EB.H for full documentation).
// ---------------------------------------------------------------
void voxelize_stl_double_precision(
    const std::string& stl_file,
    amrex::Real scale_r,
    const amrex::Array<amrex::Real, 3>& center_r,
    int reverse_normal,
    const amrex::Geometry& geom,
    amrex::iMultiFab& out,
    int comp,
    int inside_value,
    int outside_value)
{
    BL_PROFILE("LBM::voxelize_stl_double_precision()");

    const double scale = static_cast<double>(scale_r);
    const double center[3] = {
        static_cast<double>(center_r[0]), static_cast<double>(center_r[1]),
        static_cast<double>(center_r[2])};

    // ---- Step 1: root reads the STL file in double, applies transform.
    std::vector<DTri> tris;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream head(stl_file, std::ios::binary);
        if (!head.good()) {
            amrex::Abort(
                "voxelize_stl_double_precision: failed to open " + stl_file);
        }
        char hdr[6] = {0};
        head.read(hdr, 5);
        head.close();
        // Binary STLs sometimes start with "solid" too, but ASCII STLs *always*
        // do.  Follow the AMReX convention (strcmp on first 5 bytes).
        const bool is_ascii = (std::strncmp(hdr, "solid", 5) == 0);

        bool ok = is_ascii ? read_ascii_stl_double(stl_file, tris)
                           : read_binary_stl_double(stl_file, tris);
        if (!ok) {
            amrex::Abort(
                "voxelize_stl_double_precision: parse failed for " + stl_file);
        }

        // Apply v -> v * scale + center in double.
        for (auto& tri : tris) {
            for (int v = 0; v < 3; ++v) {
                tri[v][0] = tri[v][0] * scale + center[0];
                tri[v][1] = tri[v][1] * scale + center[1];
                tri[v][2] = tri[v][2] * scale + center[2];
            }
            // reverse_normal semantics from amrex::STLtools: swap v1<->v2.
            // The intersection test is winding-independent, but we honour
            // the flip for parity with the AMReX path.
            if (reverse_normal) {
                std::swap(tri[0], tri[1]);
            }
        }

        amrex::Print() << "voxelize_stl_double_precision: " << stl_file << " ("
                       << tris.size() << " triangles, scale=" << scale
                       << ", reverse_normal=" << reverse_normal << ")\n";
    }

    // ---- Step 2: broadcast the triangle list to all ranks.
    long ntri_l = static_cast<long>(tris.size());
    amrex::ParallelDescriptor::Bcast(
        &ntri_l, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    const size_t ntri = static_cast<size_t>(ntri_l);
    if (ntri == 0) {
        amrex::Abort(
            "voxelize_stl_double_precision: 0 triangles read from " + stl_file);
    }
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        tris.resize(ntri);
    }
    amrex::ParallelDescriptor::Bcast(
        reinterpret_cast<char*>(tris.data()),
        static_cast<amrex::Long>(ntri * sizeof(DTri)),
        amrex::ParallelDescriptor::IOProcessorNumber());

    // ---- Step 3: bounding box of transformed triangles + reference point.
    //
    // The ray-cast reference point must lie *outside* the mesh and must
    // not be axis-aligned with any triangle edge/vertex to avoid tie-
    // breaking headaches.  We push it well outside the bbox in a
    // deliberately non-integer direction.
    double bbmin[3] = {tris[0][0][0], tris[0][0][1], tris[0][0][2]};
    double bbmax[3] = {bbmin[0], bbmin[1], bbmin[2]};
    for (const auto& tri : tris) {
        for (int v = 0; v < 3; ++v) {
            for (int d = 0; d < 3; ++d) {
                if (tri[v][d] < bbmin[d]) bbmin[d] = tri[v][d];
                if (tri[v][d] > bbmax[d]) bbmax[d] = tri[v][d];
            }
        }
    }
    const double ptref[3] = {
        bbmin[0] - 3.11 * (bbmax[0] - bbmin[0]) - 17.13,
        bbmin[1] - 2.73 * (bbmax[1] - bbmin[1]) - 13.29,
        bbmin[2] - 2.31 * (bbmax[2] - bbmin[2]) - 7.47};

    // ---- Step 4: per-box CPU voxelization, then HtoD copy.
    const auto plo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();
    const double plo_d[3] = {
        static_cast<double>(plo[0]), static_cast<double>(plo[1]),
        static_cast<double>(plo[2])};
    const double dx_d[3] = {
        static_cast<double>(dx[0]), static_cast<double>(dx[1]),
        static_cast<double>(dx[2])};

    const amrex::IntVect nghost = out.nGrowVect();

    // Iterate over valid (owned) boxes only.  MFIter without tiling: one
    // iteration per FAB on this rank.
    for (amrex::MFIter mfi(out, false); mfi.isValid(); ++mfi) {
        const amrex::Box gbox = amrex::grow(mfi.validbox(), nghost);
        const auto lo = amrex::lbound(gbox);
        const auto hi = amrex::ubound(gbox);
        const long nx_g = static_cast<long>(hi.x - lo.x + 1);
        const long ny_g = static_cast<long>(hi.y - lo.y + 1);
        const long nz_g = static_cast<long>(hi.z - lo.z + 1);
        const long npts = nx_g * ny_g * nz_g;

        std::vector<int> hbuf(static_cast<size_t>(npts), outside_value);

#ifdef AMREX_USE_OMP
#pragma omp parallel for schedule(static)
#endif
        for (long kk = 0; kk < nz_g; ++kk) {
            const int k = static_cast<int>(lo.z + kk);
            const double zc =
                plo_d[2] + (static_cast<double>(k) + 0.5) * dx_d[2];
            for (long jj = 0; jj < ny_g; ++jj) {
                const int j = static_cast<int>(lo.y + jj);
                const double yc =
                    plo_d[1] + (static_cast<double>(j) + 0.5) * dx_d[1];
                for (long ii = 0; ii < nx_g; ++ii) {
                    const int i = static_cast<int>(lo.x + ii);
                    const double xc =
                        plo_d[0] + (static_cast<double>(i) + 0.5) * dx_d[0];

                    // Early reject if cell centre is outside the mesh bbox
                    // (with a generous margin): the segment cannot cross
                    // any triangle.  This is the same optimisation as in
                    // amrex::STLtools::fill.
                    if (xc < bbmin[0] || xc > bbmax[0] || yc < bbmin[1] ||
                        yc > bbmax[1] || zc < bbmin[2] || zc > bbmax[2]) {
                        // outside_value already initialised.
                        continue;
                    }

                    const double cellc[3] = {xc, yc, zc};
                    int crossings = 0;
                    for (const auto& tri : tris) {
                        if (segment_triangle_intersect_double(
                                ptref, cellc, tri)) {
                            ++crossings;
                        }
                    }
                    // ptref is outside -> even crossings means cell is on
                    // the same side (outside), odd means opposite (inside).
                    const int val =
                        ((crossings & 1) == 0) ? outside_value : inside_value;
                    hbuf[static_cast<size_t>(ii + nx_g * (jj + ny_g * kk))] =
                        val;
                }
            }
        }

        // Copy the host buffer into the FAB's device memory for component
        // `comp`.  BaseFab stores components contiguously in Fortran order
        // (i fastest, k slowest), so dataPtr(comp) points at the start of
        // this component's block, and the block has exactly `npts` ints.
        auto& fab = out[mfi];
        int* dst = fab.dataPtr(comp);
        amrex::Gpu::htod_memcpy(
            dst, hbuf.data(), static_cast<size_t>(npts) * sizeof(int));
    }

    amrex::Gpu::synchronize();
    out.FillBoundary(geom.periodicity());
}

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

        if constexpr (sizeof(amrex::Real) < 8) {
            // PRECISION=FLOAT: bypass amrex::STLtools (which internally
            // uses `amrex::Real` = float, and thereby produces spurious
            // solid voxels near the impeller hub in FLOAT runs, as
            // documented in EB.H::voxelize_stl_double_precision).
            // Voxelize in double directly into is_fluid.
            voxelize_stl_double_precision(
                name, scale, center, reverse_normal, geom, is_fluid,
                lbm::constants::IS_FLUID_IDX,
                /*inside_value=*/0, /*outside_value=*/1);
        } else {
            // PRECISION=DOUBLE: keep the amrex::STLtools path unchanged so
            // bit-for-bit reproducibility of prior DOUBLE plotfiles is
            // preserved.
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
                    is_fluid_arrs[nbx](i, j, k, lbm::constants::IS_FLUID_IDX) =
                        static_cast<int>(marker_arrs[nbx](i, j, k, 0));
                });
            amrex::Gpu::synchronize();
        }
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
                is_fluid_arr(i, j, k, lbm::constants::IS_FLUID_IDX) =
                    (crack_ptr[file_index] == 0) ? 1 : 0;
            });
    }
    amrex::Gpu::synchronize();

    amrex::Print() << "Voxel crack generation complete" << std::endl;
}

} // namespace lbm
