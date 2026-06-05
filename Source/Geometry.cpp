#include "Geometry.H"

namespace lbm {

void ExtrudedTriangles::build(
    const amrex::Geometry& geom, const int max_coarsening_level)
{
    // setting some constants
    // the polygon is triangle
    // we can only do a maximum of 5 triangles (change if needed)
    const int npts_in_tri = 3;
    const int max_tri = 5;

    // number of user defined triangles
    int num_tri;

    amrex::ParmParse pp("extruded_triangles");
    amrex::Vector<amrex::Array<amrex::Real, AMREX_SPACEDIM>> alltri(
        static_cast<amrex::Long>(npts_in_tri * max_tri));

    // initialize all triangles with some dummy values
    // that fall outside of the domain
    const amrex::Real* problo;
    const amrex::Real* probhi;
    amrex::Real maxlen;

    problo = geom.ProbLo();
    probhi = geom.ProbHi();

    maxlen = amrex::max<amrex::Real>(
        amrex::max<amrex::Real>(geom.ProbLength(0), geom.ProbLength(1)),
        geom.ProbLength(2));

    // setting all triangles to be waaay outside the domain initially
    for (int itri = 0; itri < max_tri; itri++) {
        alltri[npts_in_tri * itri + 0][0] = problo[0] + 100.0 * maxlen;
        alltri[npts_in_tri * itri + 0][1] = problo[1] - 100.0 * maxlen;
        alltri[npts_in_tri * itri + 0][2] = 0.0;

        alltri[npts_in_tri * itri + 1][0] = probhi[0] + 100.0 * maxlen;
        alltri[npts_in_tri * itri + 1][1] = problo[1] - 100.0 * maxlen;
        alltri[npts_in_tri * itri + 1][2] = 0.0;

        alltri[npts_in_tri * itri + 2][0] = probhi[0] + 100.0 * maxlen;
        alltri[npts_in_tri * itri + 2][1] = problo[1] + 100.0 * maxlen;
        alltri[npts_in_tri * itri + 2][2] = 0.0;
    }

    // get user defined number of triangles
    pp.get("num_tri", num_tri);

    for (int itri = 0; itri < num_tri; itri++) {
        amrex::Array<amrex::Real, AMREX_SPACEDIM> point{
            AMREX_D_DECL(0.0, 0.0, 0.0)};

        for (int ipt = 0; ipt < npts_in_tri; ipt++) {
            std::string pointstr = "tri_" + std ::to_string(itri) + "_point_" +
                                   std::to_string(ipt);
            amrex::Vector<amrex::Real> vecpt;
            pp.getarr(pointstr, vecpt, 0, AMREX_SPACEDIM);
            for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
                point[dir] = vecpt[dir];
            }
            alltri[npts_in_tri * itri + ipt] = point;
        }
    }

    // intersection of the 3 planes in a triangle for all triangles
    amrex::Vector<std::unique_ptr<amrex::EB2::IntersectionIF<
        amrex::EB2::PlaneIF, amrex::EB2::PlaneIF, amrex::EB2::PlaneIF>>>
        impfunc_triangles(max_tri);

    for (int itri = 0; itri < max_tri; itri++) {
        // make sure points are in anti clockwise direction to set the inside of
        // the triangle as solid phase correctly
        amrex::Array<amrex::Real, AMREX_SPACEDIM> norm0;
        amrex::Array<amrex::Real, AMREX_SPACEDIM> norm1;
        amrex::Array<amrex::Real, AMREX_SPACEDIM> norm2;

        amrex::Array<amrex::Real, AMREX_SPACEDIM> point0;
        amrex::Array<amrex::Real, AMREX_SPACEDIM> point1;
        amrex::Array<amrex::Real, AMREX_SPACEDIM> point2;

        point0 = alltri[npts_in_tri * itri + 0];
        point1 = alltri[npts_in_tri * itri + 1];
        point2 = alltri[npts_in_tri * itri + 2];

        norm0[0] = -(point1[1] - point0[1]);
        norm0[1] = (point1[0] - point0[0]);
        norm0[2] = 0.0;

        norm1[0] = -(point2[1] - point1[1]);
        norm1[1] = (point2[0] - point1[0]);
        norm1[2] = 0.0;

        norm2[0] = -(point0[1] - point2[1]);
        norm2[1] = (point0[0] - point2[0]);
        norm2[2] = 0.0;

        // normalize so that magnitude is 1
        amrex::Real norm = sqrt(norm0[0] * norm0[0] + norm0[1] * norm0[1]);
        norm0[0] = norm0[0] / norm;
        norm0[1] = norm0[1] / norm;

        // normalize so that magnitude is 1
        norm = sqrt(norm1[0] * norm1[0] + norm1[1] * norm1[1]);
        norm1[0] = norm1[0] / norm;
        norm1[1] = norm1[1] / norm;

        // normalize so that magnitude is 1
        norm = sqrt(norm2[0] * norm2[0] + norm2[1] * norm2[1]);
        norm2[0] = norm2[0] / norm;
        norm2[1] = norm2[1] / norm;

        amrex::EB2::PlaneIF plane0(point0, norm0);
        amrex::EB2::PlaneIF plane1(point1, norm1);
        amrex::EB2::PlaneIF plane2(point2, norm2);

        impfunc_triangles[itri] = std::make_unique<amrex::EB2::IntersectionIF<
            amrex::EB2::PlaneIF, amrex::EB2::PlaneIF, amrex::EB2::PlaneIF>>(

            plane0, plane1, plane2);
    }

    auto alltri_if = amrex::EB2::makeUnion(
        *impfunc_triangles[0], *impfunc_triangles[1], *impfunc_triangles[2],
        *impfunc_triangles[3], *impfunc_triangles[4]);

    auto alltri_extrude_if = amrex::EB2::extrude(alltri_if, 2); // along z

    auto gshop = amrex::EB2::makeShop(alltri_extrude_if);
    amrex::EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level);
}

void RotatedCylinder::build(
    const amrex::Geometry& geom, const int max_coarsening_level)
{
    bool inside = true;
    amrex::Real radius = 0.0002;
    int direction = 0;
    amrex::Vector<amrex::Real> centervec(3, 0.0);
    amrex::Real rotation = 0;
    int rotation_axe = 0;

    amrex::ParmParse pp("eb2");
    pp.query("cylinder_has_fluid_inside", inside);
    pp.query("cylinder_radius", radius);
    pp.query("cylinder_direction", direction);
    pp.query("cylinder_rotation", rotation);
    pp.query("cylinder_rotation_axe", rotation_axe);
    pp.getarr("cylinder_center", centervec, 0, 3);
    amrex::Array<amrex::Real, AMREX_SPACEDIM> center = {
        AMREX_D_DECL(centervec[0], centervec[1], centervec[2])};

    rotation = (rotation / 180.) * M_PI;

    amrex::EB2::CylinderIF my_cyl(radius, direction, center, inside);

    auto my_cyl_rot = amrex::EB2::rotate(my_cyl, rotation, rotation_axe);
    auto gshop = amrex::EB2::makeShop(my_cyl_rot);
    amrex::EB2::Build(
        gshop, geom, max_coarsening_level, max_coarsening_level, 4, false);
}

void RotatedBox::build(
    const amrex::Geometry& geom, const int max_coarsening_level)
{

    amrex::ParmParse pp("eb2");
    amrex::RealArray lo;
    pp.get("box_lo", lo);

    amrex::RealArray hi;
    pp.get("box_hi", hi);

    bool has_fluid_inside;
    pp.get("box_has_fluid_inside", has_fluid_inside);

    amrex::Real rotation = 0;
    int rotation_axe = 0;
    amrex::RealArray rotation_center({0.0});
    pp.query("box_rotation", rotation);
    pp.query("box_rotation_axe", rotation_axe);
    pp.query("box_rotation_center", rotation_center);
    const amrex::RealArray neg_rotation_center({AMREX_D_DECL(
        -rotation_center[0], -rotation_center[1], -rotation_center[2])});

    rotation = (rotation / 180.) * M_PI;

    amrex::EB2::BoxIF bf(lo, hi, has_fluid_inside);

    auto bf_rot = amrex::EB2::translate(
        amrex::EB2::rotate(
            amrex::EB2::translate(bf, neg_rotation_center), rotation,
            rotation_axe),
        rotation_center);
    auto gshop = amrex::EB2::makeShop(bf_rot);
    amrex::EB2::Build(
        gshop, geom, max_coarsening_level, max_coarsening_level, 4, false);
}

} // namespace lbm
