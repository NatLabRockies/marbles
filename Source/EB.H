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

} // namespace lbm
#endif
