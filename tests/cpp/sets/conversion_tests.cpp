#include "../utils.h"
#include "catch2/catch_test_macros.hpp"
#include <Eigen/Dense>
#include <geosets/conversion/to_hpolytope.h>
#include <geosets/conversion/to_polygon.h>
#include <geosets/conversion/to_vpolytope.h>
#include <geosets/sets/polygon.h>
#include <geosets/sets/vpolytope.h>
#include <geosets/sets/zonotope.h>

using namespace geosets;

namespace {

constexpr double tol = 1e-8;

} // namespace

TEST_CASE("VPolytope to HPolytope conversion in full dimension", "[conversion][to_hpolytope]") {
    Eigen::MatrixXd vertices(4, 2);
    vertices << -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0;

    const VPolytope poly(vertices);
    const HPolytope hpoly = conversion::vpolytope_to_hpolytope(poly);

    REQUIRE(hpoly.dim() == 2);
    REQUIRE(hpoly.n_constraints() >= 4);
    REQUIRE(hpoly.contains_point(Eigen::Vector2d::Zero(), tol));
    REQUIRE(hpoly.contains_point(Eigen::Vector2d(-1.0, 1.0), tol));
    REQUIRE_FALSE(hpoly.contains_point(Eigen::Vector2d(1.5, 0.0), tol));
}

TEST_CASE("VPolytope to HPolytope conversion for degenerate line in 2D",
          "[conversion][to_hpolytope][degenerate]") {
    Eigen::MatrixXd vertices(3, 2);
    vertices << -2.0, -1.0, 0.0, 0.0, 2.0, 1.0;

    const VPolytope poly(vertices, false);
    const HPolytope hpoly = conversion::vpolytope_to_hpolytope(poly);

    REQUIRE(hpoly.dim() == 2);
    REQUIRE(hpoly.contains_point(Eigen::Vector2d(-2.0, -1.0), tol));
    REQUIRE(hpoly.contains_point(Eigen::Vector2d(0.0, 0.0), tol));
    REQUIRE(hpoly.contains_point(Eigen::Vector2d(2.0, 1.0), tol));
    REQUIRE_FALSE(hpoly.contains_point(Eigen::Vector2d(0.0, 0.2), tol));
    REQUIRE_FALSE(hpoly.contains_point(Eigen::Vector2d(-3.0, -1.5), tol));
}

TEST_CASE("VPolytope to HPolytope conversion for single vertex",
          "[conversion][to_hpolytope][vertex]") {
    Eigen::MatrixXd vertices(1, 3);
    vertices << -1.0, 2.0, -3.0;

    const VPolytope poly(vertices, false);
    const HPolytope hpoly = conversion::vpolytope_to_hpolytope(poly);

    REQUIRE(hpoly.dim() == 3);
    REQUIRE(hpoly.contains_point(Eigen::Vector3d(-1.0, 2.0, -3.0), tol));
    REQUIRE_FALSE(hpoly.contains_point(Eigen::Vector3d(-1.0, 2.0, -2.9), tol));
}

TEST_CASE("VPolytope to HPolytope conversion for planar set in 3D",
          "[conversion][to_hpolytope][affine]") {
    Eigen::MatrixXd vertices(4, 3);
    vertices << -2.0, -1.0, -2.0, 2.0, -1.0, -2.0, 2.0, 1.0, -2.0, -2.0, 1.0, -2.0;

    const VPolytope poly(vertices);
    const HPolytope hpoly = conversion::vpolytope_to_hpolytope(poly);

    REQUIRE(hpoly.dim() == 3);
    REQUIRE(hpoly.contains_point(Eigen::Vector3d(0.0, 0.0, -2.0), tol));
    REQUIRE(hpoly.contains_point(Eigen::Vector3d(-2.0, 1.0, -2.0), tol));
    REQUIRE_FALSE(hpoly.contains_point(Eigen::Vector3d(0.0, 0.0, -1.0), tol));
    REQUIRE_FALSE(hpoly.contains_point(Eigen::Vector3d(3.0, 0.0, -2.0), tol));
}

TEST_CASE("HPolytope to VPolytope conversion for unit box", "[conversion][to_vpolytope]") {
    const HPolytope hpoly = HPolytope::from_unit_box(2);
    const VPolytope vpoly = conversion::hpolytope_to_vpolytope(hpoly);

    Eigen::MatrixXd expected(4, 2);
    expected << -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0;

    REQUIRE(vertices_to_set(vpoly.vertices) == vertices_to_set(expected));
}

TEST_CASE("HPolytope to VPolytope conversion for degenerate segment",
          "[conversion][to_vpolytope][degenerate]") {
    Eigen::MatrixXd A(6, 2);
    Eigen::VectorXd b(6);
    A << 1.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 0.0, -1.0, 0.0;
    b << 2.0, 2.0, 0.0, 0.0, 2.0, 2.0;

    const HPolytope hpoly(A, b);
    const VPolytope vpoly = conversion::hpolytope_to_vpolytope(hpoly);

    Eigen::MatrixXd expected(2, 2);
    expected << -2.0, 0.0, 2.0, 0.0;

    REQUIRE(vertices_to_set(vpoly.vertices) == vertices_to_set(expected));
}

TEST_CASE("HPolytope::to_vertices uses H-to-V conversion", "[conversion][to_vertices]") {
    const HPolytope hpoly = HPolytope::from_unit_box(2);
    const Eigen::MatrixXd vertices = hpoly.to_vertices();

    Eigen::MatrixXd expected(4, 2);
    expected << -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0;

    REQUIRE(vertices_to_set(vertices) == vertices_to_set(expected));
}

TEST_CASE("Polygon to VPolytope conversion for unit box", "[conversion][polygon_to_vpolytope]") {
    const Polygon polygon = Polygon::from_unit_box(2);
    const VPolytope vpoly = conversion::polygon_to_vpolytope(polygon);

    Eigen::MatrixXd expected(4, 2);
    expected << -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0;

    REQUIRE(vertices_to_set(vpoly.vertices) == vertices_to_set(expected));
}

TEST_CASE("VPolytope to Polygon conversion uses convex hull",
          "[conversion][vpolytope_to_polygon]") {
    Eigen::MatrixXd vertices(6, 2);
    vertices << 1.0, -1.0, -1.0, -1.0, 0.0, 0.0, -1.0, 1.0, 1.0, 1.0, 0.0, 0.2;

    const VPolytope vpoly(vertices);
    const Polygon polygon = conversion::vpolytope_to_polygon(vpoly);

    Eigen::MatrixXd expected(4, 2);
    expected << -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0;

    REQUIRE(vertices_to_set(polygon.to_vertices()) == vertices_to_set(expected));
}

TEST_CASE("VPolytope to Polygon conversion rejects non-2D sets",
          "[conversion][vpolytope_to_polygon][dimension]") {
    const VPolytope box3d = VPolytope::from_unit_box(3);
    REQUIRE_THROWS_AS(conversion::vpolytope_to_polygon(box3d), PolygonDimensionException);
}

// --- Zonotope to HPolytope ---

// Helper: check that all zonotope vertices are inside the HPolytope,
// and that points slightly outside (scaled from center) are not.
static void check_zonotope_hpolytope_consistency(const Zonotope& zono, const HPolytope& hpoly,
                                                 double tol_in = tol, double overshoot = 1.5) {
    REQUIRE(hpoly.dim() == zono.dim());
    Eigen::VectorXd c = zono.center();

    // All vertices must be inside
    Eigen::MatrixXd verts = zono.to_vertices();
    for (int i = 0; i < verts.rows(); ++i) {
        REQUIRE(hpoly.contains_point(verts.row(i).transpose(), tol_in));
    }

    // Overshoot each vertex from center — should be outside
    for (int i = 0; i < verts.rows(); ++i) {
        Eigen::VectorXd outside = c + overshoot * (verts.row(i).transpose() - c);
        REQUIRE_FALSE(hpoly.contains_point(outside, tol_in));
    }
}

TEST_CASE("Zonotope to HPolytope full-rank", "[conversion][zonotope_to_hpolytope]") {
    SECTION("unit box 2D") {
        check_zonotope_hpolytope_consistency(
            Zonotope::from_unit_box(2),
            conversion::zonotope_to_hpolytope(Zonotope::from_unit_box(2)));
    }
    SECTION("unit box 3D") {
        check_zonotope_hpolytope_consistency(
            Zonotope::from_unit_box(3),
            conversion::zonotope_to_hpolytope(Zonotope::from_unit_box(3)));
    }
    SECTION("offset center") {
        Eigen::Vector2d c(2.0, 3.0);
        Eigen::MatrixXd G(2, 2);
        G << 1.0, 0.0, 0.0, 1.0;
        Zonotope zono(c, G);
        check_zonotope_hpolytope_consistency(zono, conversion::zonotope_to_hpolytope(zono));
    }
    SECTION("more generators than dimensions") {
        Eigen::Vector2d c(0.0, 0.0);
        Eigen::MatrixXd G(2, 3);
        G << 1.0, 0.0, 1.0, 0.0, 1.0, 1.0;
        Zonotope zono(c, G);
        check_zonotope_hpolytope_consistency(zono, conversion::zonotope_to_hpolytope(zono));
    }
}

TEST_CASE("Zonotope to HPolytope degenerate", "[conversion][zonotope_to_hpolytope][degenerate]") {
    SECTION("line in 2D") {
        Eigen::Vector2d c(1.0, 2.0);
        Eigen::MatrixXd G(2, 1);
        G << 1.0, 0.0;
        Zonotope zono(c, G);
        check_zonotope_hpolytope_consistency(zono, conversion::zonotope_to_hpolytope(zono));
    }
    SECTION("plane in 3D") {
        Eigen::Vector3d c(0.0, 0.0, 1.0);
        Eigen::MatrixXd G(3, 2);
        G << 1.0, 0.0, 0.0, 1.0, 0.0, 0.0;
        Zonotope zono(c, G);
        check_zonotope_hpolytope_consistency(zono, conversion::zonotope_to_hpolytope(zono));
    }
    SECTION("single point") {
        Eigen::Vector3d c(1.0, 2.0, 3.0);
        Eigen::MatrixXd G = Eigen::MatrixXd::Zero(3, 2);
        Zonotope zono(c, G);
        HPolytope hpoly = conversion::zonotope_to_hpolytope(zono);
        REQUIRE(hpoly.dim() == 3);
        REQUIRE(hpoly.contains_point(c, tol));
        REQUIRE_FALSE(hpoly.contains_point(c + Eigen::Vector3d(0.1, 0.0, 0.0), tol));
    }
}
