#ifndef EB_H
#define EB_H
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Scale.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB2_IF_Rotation.H>
#include <AMReX_EB2_IF_Lathe.H>
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Ellipsoid.H>
#include <AMReX_EB2_IF_Sphere.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Extrusion.H>
#include <AMReX_EB2_GeometryShop.H>
#include <AMReX_EB2_IF_Difference.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBSupport.H>
#include <AMReX_EB_STL_utils.H>

namespace lbm {
void initialize_eb(const amrex::Geometry& geom, const int max_level);
void initialize_from_stl(
    const amrex::Geometry& geom, amrex::iMultiFab& is_fluid);
std::vector<uint16_t>
read_crack_file(const std::string& filename, int nx, int ny, int nz);

void generate_voxel_cracks(
    const amrex::Geometry& geom, amrex::iMultiFab& is_fluid);

// Fill component `comp` of `out` with `outside_value` (default 1) for
// cells whose center lies outside the closed surface described by `stl_file`,
// and with `inside_value` (default 0) for cells whose center lies inside.
// The STL vertices are transformed by v -> v * scale + center before
// classification.  `reverse_normal` != 0 flips the mesh winding (v2 <-> v3
// after transformation), matching amrex::STLtools::read_stl_file semantics.
//
// The rasterization runs entirely in double precision on the host regardless
// of the compile-time PRECISION of `amrex::Real`, then the resulting int
// mask is HtoD-copied into `out`.  This exists because amrex::STLtools uses
// `amrex::Real` internally; in a PRECISION=FLOAT build the ray-triangle
// tests are done in `float` (eps ~= 1.2e-7), which produces spurious
// classification flips at cells whose center lies within eps of the closed
// surface -- most visibly near the impeller hub where the STL has many
// small triangles meeting at the shaft.  This routine reads the STL as
// `float32` (STL binary format is 4-byte floats by definition), widens to
// double for every arithmetic op, and thereby matches the DOUBLE-build
// result cell-for-cell (verified against resultsJuly13/plt00000 for
// pitched_impeller.stl on the 180^3 kLa bioreactor grid).
//
// Ghost cells covered by the box array's periodicity are also written.
// If `stl_file` cannot be opened or contains no triangles, aborts.
void voxelize_stl_double_precision(
    const std::string& stl_file, amrex::Real scale,
    const amrex::Array<amrex::Real, 3>& center, int reverse_normal,
    const amrex::Geometry& geom, amrex::iMultiFab& out, int comp = 0,
    int inside_value = 0, int outside_value = 1);

} // namespace lbm
#endif
