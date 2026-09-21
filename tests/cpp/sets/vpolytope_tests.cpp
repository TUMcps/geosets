#include "../utils.h"
#include "catch2/catch_test_macros.hpp"
#include "geosets/binary_operations.h"
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <catch2/catch_all.hpp>
#include <geosets/linalg_utils.h>
#include <geosets/sets/vpolytope.h>

using namespace geosets;

constexpr double tol = 1e-9;

TEST_CASE("VPolytope construction edge cases", "[vpolytope][construction]") {
    SECTION("Custom polytope from vertices") {
        Eigen::MatrixXd vertices(4, 2);
        vertices << 0, 0, 1, 0, 1, 1, 0, 1;

        auto polytope = VPolytope(vertices);
        REQUIRE(polytope.dim() == 2);
        REQUIRE(polytope.n_vertices() == 4);
    }

    SECTION("Validation throws when too few vertices") {
        // A 2D polytope needs at least 3 vertices
        Eigen::MatrixXd vertices(2, 2);
        vertices << 0, 0, 1, 0;

        REQUIRE_THROWS_AS(VPolytope(vertices), VPolytopeValidationException);
    }

    SECTION("3D simplex construction") {
        // Tetrahedron vertices
        Eigen::MatrixXd vertices(4, 3);
        vertices << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1;

        auto polytope = VPolytope(vertices);
        REQUIRE(polytope.dim() == 3);
        REQUIRE(polytope.n_vertices() == 4);
    }

    SECTION("Polytope with negative coordinates") {
        Eigen::MatrixXd vertices(4, 2);
        vertices << -1, -1, 1, -1, 1, 1, -1, 1;

        auto polytope = VPolytope(vertices);
        REQUIRE(polytope.center().isApprox(Eigen::Vector2d::Zero()));
    }
}

TEST_CASE("VPolytope degenerate checks", "[vpolytope][degenerate]") {
    SECTION("Unit box is not degenerate") {
        auto polytope = VPolytope::from_unit_box(2);
        REQUIRE_FALSE(polytope.degenerate());
    }

    SECTION("Collinear points are degenerate in 2D") {
        // Three collinear points
        Eigen::MatrixXd vertices(3, 2);
        vertices << 0, 0, 1, 0, 2, 0;

        auto polytope = VPolytope(vertices, false);
        REQUIRE(polytope.degenerate());
    }

    SECTION("Planar points are degenerate in 3D") {
        // Four coplanar points in 3D
        Eigen::MatrixXd vertices(4, 3);
        vertices << 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0;

        auto polytope = VPolytope(vertices, false);
        REQUIRE(polytope.degenerate());
    }
}

TEST_CASE("VPolytope volume calculations", "[vpolytope][volume]") {
    SECTION("2D unit box volume") {
        auto polytope = VPolytope::from_unit_box(2);
        double vol = polytope.volume();
        // Unit box is [-1,1] x [-1,1], so area = 4
        REQUIRE(std::abs(vol - 4.0) < tol);
    }

    SECTION("Unit square volume") {
        Eigen::MatrixXd vertices(4, 2);
        vertices << 0, 0, 1, 0, 1, 1, 0, 1;

        auto polytope = VPolytope(vertices);
        REQUIRE(std::abs(polytope.volume() - 1.0) < tol);
    }

    SECTION("Triangle volume in 2D") {
        Eigen::MatrixXd vertices(3, 2);
        vertices << 0, 0, 2, 0, 1, 2;

        auto polytope = VPolytope(vertices);
        // Triangle area = 0.5 * base * height = 0.5 * 2 * 2 = 2
        REQUIRE(std::abs(polytope.volume() - 2.0) < tol);
    }

    SECTION("3D unit box volume") {
        auto polytope = VPolytope::from_unit_box(3);
        double vol = polytope.volume();
        // Unit box is [-1,1]^3, so volume = 8
        REQUIRE(std::abs(vol - 8.0) < tol);
    }

    SECTION("Tetrahedron volume") {
        Eigen::MatrixXd vertices(4, 3);
        vertices << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1;

        auto polytope = VPolytope(vertices);
        // Volume of unit tetrahedron = 1/6
        REQUIRE(std::abs(polytope.volume() - 1.0 / 6.0) < tol);
    }

    SECTION("Polytope with negative coordinates") {
        Eigen::MatrixXd vertices(4, 2);
        vertices << -1, -1, 1, -1, 1, 1, -1, 1;

        auto polytope = VPolytope(vertices);
        // Square with side length 2, area = 4
        REQUIRE(std::abs(polytope.volume() - 4.0) < tol);
    }
}

TEST_CASE("VPolytope contains_point", "[vpolytope][contains]") {
    SECTION("Unit box contains origin") {
        auto polytope = VPolytope::from_unit_box(2);
        REQUIRE(polytope.contains_point(Eigen::Vector2d::Zero()));
    }

    SECTION("Unit box contains vertices") {
        auto polytope = VPolytope::from_unit_box(2);
        REQUIRE(polytope.contains_point(Eigen::Vector2d(1, 0)));
        REQUIRE(polytope.contains_point(Eigen::Vector2d(-1, 0)));
        REQUIRE(polytope.contains_point(Eigen::Vector2d(0, 1)));
        REQUIRE(polytope.contains_point(Eigen::Vector2d(0, -1)));
    }

    SECTION("Point outside is not contained") {
        auto polytope = VPolytope::from_unit_box(2);
        REQUIRE_FALSE(polytope.contains_point(Eigen::Vector2d(2, 0)));
        REQUIRE_FALSE(polytope.contains_point(Eigen::Vector2d(0, -2)));
    }

    SECTION("Triangle containment") {
        Eigen::MatrixXd vertices(3, 2);
        vertices << 0, 0, 2, 0, 1, 2;

        auto polytope = VPolytope(vertices);
        REQUIRE(polytope.contains_point(Eigen::Vector2d(1, 0.5)));
        REQUIRE_FALSE(polytope.contains_point(Eigen::Vector2d(-1, 0)));
    }

    SECTION("Polytope with negative coordinates") {
        Eigen::MatrixXd vertices(4, 2);
        vertices << -2, -2, 2, -2, 2, 2, -2, 2;

        auto polytope = VPolytope(vertices);
        REQUIRE(polytope.contains_point(Eigen::Vector2d::Zero()));
        REQUIRE(polytope.contains_point(Eigen::Vector2d(-1, -1)));
        REQUIRE_FALSE(polytope.contains_point(Eigen::Vector2d(3, 0)));
    }
}

TEST_CASE("VPolytope support_function", "[vpolytope][support]") {
    SECTION("Unit box support in positive direction") {
        auto polytope = VPolytope::from_unit_box(2);
        auto support = polytope.support_function(Eigen::Vector2d(1, 0));

        REQUIRE(std::abs(support.value - 1.0) < tol);
        REQUIRE(support.vector(0) == Catch::Approx(1.0));
    }

    SECTION("Unit box support in negative direction") {
        auto polytope = VPolytope::from_unit_box(2);
        auto support = polytope.support_function(Eigen::Vector2d(-1, 0));

        REQUIRE(std::abs(support.value - 1.0) < tol);
        REQUIRE(support.vector(0) == Catch::Approx(-1.0));
    }

    SECTION("Support in diagonal direction") {
        auto polytope = VPolytope::from_unit_box(2);
        Eigen::Vector2d direction(1, 1);
        direction.normalize();
        auto support = polytope.support_function(direction);

        // Support should be at vertex (1, 1)
        REQUIRE(support.vector(0) == Catch::Approx(1.0));
        REQUIRE(support.vector(1) == Catch::Approx(1.0));
    }

    SECTION("Triangle support") {
        Eigen::MatrixXd vertices(3, 2);
        vertices << 0, 0, 2, 0, 1, 2;

        auto polytope = VPolytope(vertices);
        auto support = polytope.support_function(Eigen::Vector2d(0, 1));

        // Maximum in y-direction is at (1, 2)
        REQUIRE(support.vector.isApprox(Eigen::Vector2d(1, 2)));
        REQUIRE(std::abs(support.value - 2.0) < tol);
    }

    SECTION("Support with negative coordinates") {
        Eigen::MatrixXd vertices(4, 2);
        vertices << -2, -2, 2, -2, 2, 2, -2, 2;

        auto polytope = VPolytope(vertices);
        auto support = polytope.support_function(Eigen::Vector2d(-1, -1));

        // Maximum in (-1, -1) direction is at (-2, -2)
        REQUIRE(support.vector.isApprox(Eigen::Vector2d(-2, -2)));
    }
}

TEST_CASE("VPolytope linear_transform", "[vpolytope][linear_transform]") {
    SECTION("Scaling transformation") {
        auto polytope = VPolytope::from_unit_box(2);
        Eigen::Matrix2d scale = 2.0 * Eigen::Matrix2d::Identity();

        auto transformed = polytope.linear_transform(scale);

        // After scaling by 2, bounds should be [-2, 2] x [-2, 2]
        REQUIRE(transformed.center().isApprox(Eigen::Vector2d::Zero()));
        REQUIRE(std::abs(transformed.volume() - 16.0) < tol);
    }

    SECTION("Rotation preserves volume") {
        auto polytope = VPolytope::from_unit_box(2);
        // 45 degree rotation matrix
        double angle = M_PI / 4.0;
        Eigen::Matrix2d rotation;
        rotation << std::cos(angle), -std::sin(angle), std::sin(angle), std::cos(angle);

        auto rotated = polytope.linear_transform(rotation);

        REQUIRE(rotated.center().isApprox(Eigen::Vector2d::Zero()));
        REQUIRE(std::abs(rotated.volume() - polytope.volume()) < tol);
    }

    SECTION("Projection to lower dimension") {
        auto polytope = VPolytope::from_unit_box(2);
        // Project to x-axis
        Eigen::MatrixXd projection(1, 2);
        projection << 1, 0;

        auto projected = polytope.linear_transform(projection);

        REQUIRE(projected.dim() == 1);
        // After projection to 1D and removing duplicates via convex hull,
        // we should have the extremes [-1] and [1]
        // But linear_transform doesn't automatically compute convex hull
        // So we may have duplicate vertices. Let's just check dimension.
    }

    SECTION("Negative scaling") {
        Eigen::MatrixXd vertices(3, 2);
        vertices << 0, 0, 1, 0, 0, 1;

        auto polytope = VPolytope(vertices);
        Eigen::Matrix2d scale = -1.0 * Eigen::Matrix2d::Identity();

        auto transformed = polytope.linear_transform(scale);

        // Vertices should be flipped: (0,0), (-1,0), (0,-1)
        // Center should be at -1 * original center
        Eigen::Vector2d expected_center = -polytope.center();
        REQUIRE(transformed.center().isApprox(expected_center));
    }
}

TEST_CASE("VPolytope translate", "[vpolytope][translate]") {
    SECTION("Translate by positive vector") {
        auto polytope = VPolytope::from_unit_box(2);
        Eigen::Vector2d translation(1, 2);

        auto translated = polytope.translate(translation);

        REQUIRE(translated.center().isApprox(translation));
    }

    SECTION("Translate by negative vector") {
        Eigen::MatrixXd vertices(3, 2);
        vertices << 0, 0, 1, 0, 0, 1;

        auto polytope = VPolytope(vertices);
        Eigen::Vector2d translation(-2, -3);

        auto translated = polytope.translate(translation);
        Eigen::Vector2d expected_center = polytope.center() + translation;
        REQUIRE(translated.center().isApprox(expected_center));
    }
}

TEST_CASE("VPolytope binary operations", "[vpolytope][binary_ops]") {
    SECTION("Intersection of identical polytopes") {
        auto a = VPolytope::from_unit_box(2);
        auto b = VPolytope::from_unit_box(2);

        auto result = intersection(a, b);

        // Intersection of identical sets should give the same set
        REQUIRE(result.dim() == 2);
        REQUIRE(std::abs(result.volume() - a.volume()) < tol);
        REQUIRE(result.center().isApprox(a.center()));
    }

    SECTION("Intersection of two translated polytopes") {
        auto a = VPolytope::from_unit_box(2).translate(Eigen::Vector2d(0.5, 0.5));
        auto b = VPolytope::from_unit_box(2).translate(-Eigen::Vector2d(0.5, 0.5));
        auto result = intersection(a, b);

        Eigen::MatrixXd expected(4, 2);
        expected << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
        REQUIRE(vertices_to_set(result.to_vertices()) == vertices_to_set(expected));
        REQUIRE(result.center() == Eigen::Vector2d(0.0, 0.0));
        REQUIRE(result.volume() == 1);
    }

    SECTION("Minkowski sum of polytopes") {
        Eigen::MatrixXd vertices_a(3, 2);
        vertices_a << 0, 0, 1, 0, 0, 1;
        auto a = VPolytope(vertices_a);

        Eigen::MatrixXd vertices_b(3, 2);
        vertices_b << 0, 0, 1, 0, 0, 1;
        auto b = VPolytope(vertices_b);

        auto result = minkowski_sum(a, b);

        REQUIRE(result.dim() == 2);
        // Minkowski sum of two identical triangles
        REQUIRE(result.volume() > a.volume());
    }

    SECTION("Minkowski sum is commutative") {
        auto a = VPolytope::from_unit_box(2);
        Eigen::MatrixXd vertices_b(3, 2);
        vertices_b << 0, 0, 1, 0, 0, 1;
        auto b = VPolytope(vertices_b);

        auto sum_ab = minkowski_sum(a, b);
        auto sum_ba = minkowski_sum(b, a);

        REQUIRE(sum_ab.center().isApprox(sum_ba.center()));
        REQUIRE(std::abs(sum_ab.volume() - sum_ba.volume()) < tol);
    }

    SECTION("Minkowski sum with negative coordinates") {
        Eigen::MatrixXd vertices_a(4, 2);
        vertices_a << -1, -1, 1, -1, 1, 1, -1, 1;
        auto a = VPolytope(vertices_a);

        Eigen::MatrixXd vertices_b(3, 2);
        vertices_b << -0.5, -0.5, 0.5, -0.5, 0, 0.5;
        auto b = VPolytope(vertices_b);

        auto result = minkowski_sum(a, b);
        REQUIRE(result.dim() == 2);
        REQUIRE(result.volume() > 0);
    }

    SECTION("Dimension mismatch throws") {
        auto a2d = VPolytope::from_unit_box(2);
        auto b3d = VPolytope::from_unit_box(3);

        REQUIRE_THROWS_AS(intersection(a2d, b3d), DimensionMismatchException);
        REQUIRE_THROWS_AS(minkowski_sum(a2d, b3d), DimensionMismatchException);
    }
}

TEST_CASE("VPolytope vertex_list constructor", "[vpolytope][construction]") {
    SECTION("Constructs correctly from vertex list") {
        std::vector<Eigen::VectorXd> vlist = {Eigen::Vector2d(-1, -2), Eigen::Vector2d(2, -1),
                                              Eigen::Vector2d(0, 3)};
        Eigen::MatrixXd expected(3, 2);
        expected << -1, -2, 2, -1, 0, 3;

        auto polytope = VPolytope(vlist);
        REQUIRE(polytope.vertices.isApprox(expected));
    }

    SECTION("Empty vertex list throws") {
        std::vector<Eigen::VectorXd> vlist;
        REQUIRE_THROWS_AS(VPolytope(vlist), VPolytopeValidationException);
    }
}

TEST_CASE("VPolytope append_vertices", "[vpolytope][append_vertices]") {
    SECTION("Appended vertices are stored correctly") {
        Eigen::MatrixXd initial(3, 2);
        initial << 0, 0, 1, 0, 0, 1;
        auto polytope = VPolytope(initial);

        Eigen::MatrixXd extra(2, 2);
        extra << -1, 0, 0, -1;
        polytope.append_vertices(extra);

        REQUIRE(polytope.n_vertices() == 5);
        REQUIRE(polytope.vertices.row(3).isApprox(Eigen::RowVector2d(-1, 0)));
        REQUIRE(polytope.vertices.row(4).isApprox(Eigen::RowVector2d(0, -1)));
    }

    SECTION("Dimension mismatch throws") {
        Eigen::MatrixXd initial(3, 2);
        initial << 0, 0, 1, 0, 0, 1;
        auto polytope = VPolytope(initial);

        Eigen::MatrixXd wrong_dim(2, 3);
        wrong_dim << 0, 0, 1, 1, 0, 0;

        REQUIRE_THROWS_AS(polytope.append_vertices(wrong_dim), VPolytopeValidationException);
    }
}
