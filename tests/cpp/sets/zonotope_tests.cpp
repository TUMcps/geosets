#include "catch2/catch_test_macros.hpp"
#include <Eigen/Dense>
#include <catch2/catch_all.hpp>
#include <cmath>
#include <geosets/binary_operations.h>
#include <geosets/sets/zonotope.h>

using namespace geosets;

TEST_CASE("Zonotope construction", "[zonotope]") {
    Eigen::VectorXd c(2);
    c << 1.0, 2.0;
    Eigen::MatrixXd G(2, 3);
    G << 1, 0, 1, 0, 1, 1;

    Zonotope z(c, G);
    REQUIRE(z.dim() == 2);
    REQUIRE(z.n_generators() == 3);
    REQUIRE(z.order() == Catch::Approx(1.5));
}

TEST_CASE("Zonotope validation", "[zonotope]") {
    Eigen::VectorXd c(3);
    c << 1, 2, 3;
    Eigen::MatrixXd G(2, 2);
    G << 1, 0, 0, 1;

    REQUIRE_THROWS_AS(Zonotope(c, G), ZonotopeValidationException);
}

TEST_CASE("Zonotope center", "[zonotope]") {
    Eigen::VectorXd c(2);
    c << 3.0, 4.0;
    Eigen::MatrixXd G = Eigen::MatrixXd::Identity(2, 2);

    Zonotope z(c, G);
    REQUIRE(z.center().isApprox(c));
}

TEST_CASE("Zonotope support function", "[zonotope]") {
    // Unit box zonotope: c=0, G=I
    auto z = Zonotope::from_unit_box(2);

    SECTION("axis-aligned direction") {
        Eigen::VectorXd d(2);
        d << 1.0, 0.0;
        auto sup = z.support_function(d);
        REQUIRE(sup.value == Catch::Approx(1.0));
        REQUIRE(sup.vector.isApprox(Eigen::Vector2d(1.0, 1.0), 1e-9));
    }

    SECTION("diagonal direction") {
        Eigen::VectorXd d(2);
        d << 1.0, 1.0;
        auto sup = z.support_function(d);
        REQUIRE(sup.value == Catch::Approx(2.0));
    }
}

TEST_CASE("Zonotope contains point", "[zonotope]") {
    auto z = Zonotope::from_unit_box(2);

    REQUIRE(z.contains_point(Eigen::Vector2d(0.0, 0.0)));
    REQUIRE(z.contains_point(Eigen::Vector2d(1.0, 1.0)));
    REQUIRE(z.contains_point(Eigen::Vector2d(-1.0, -1.0)));
    REQUIRE(z.contains_point(Eigen::Vector2d(0.5, -0.3)));
    REQUIRE_FALSE(z.contains_point(Eigen::Vector2d(1.5, 0.0)));
    REQUIRE_FALSE(z.contains_point(Eigen::Vector2d(0.0, 1.5)));
}

TEST_CASE("Zonotope contains point - non-axis-aligned", "[zonotope]") {
    Eigen::VectorXd c(2);
    c << 0.0, 0.0;
    Eigen::MatrixXd G(2, 2);
    G << 1, 1, 0, 1;

    Zonotope z(c, G);

    // c + G*[1,1] = [2,1], c + G*[-1,1] = [0,1], etc.
    REQUIRE(z.contains_point(Eigen::Vector2d(0.0, 0.0)));
    REQUIRE(z.contains_point(Eigen::Vector2d(2.0, 1.0)));
    REQUIRE_FALSE(z.contains_point(Eigen::Vector2d(3.0, 0.0)));
}

TEST_CASE("Zonotope volume - unit box", "[zonotope]") {
    SECTION("2D") {
        auto z = Zonotope::from_unit_box(2);
        REQUIRE(z.volume() == Catch::Approx(4.0)); // 2^2 * det(I) = 4
    }
    SECTION("3D") {
        auto z = Zonotope::from_unit_box(3);
        REQUIRE(z.volume() == Catch::Approx(8.0)); // 2^3 * det(I) = 8
    }
}

TEST_CASE("Zonotope volume - with extra generators", "[zonotope]") {
    // 2D zonotope with 3 generators: hexagon-like
    Eigen::VectorXd c = Eigen::Vector2d::Zero();
    Eigen::MatrixXd G(2, 3);
    G << 1, 0, 1, 0, 1, 1;

    Zonotope z(c, G);
    // V = 2^2 * (|det([g1,g2])| + |det([g1,g3])| + |det([g2,g3])|)
    // = 4 * (|1| + |1| + |-1|) = 12
    REQUIRE(z.volume() == Catch::Approx(12.0));
}

TEST_CASE("Zonotope degenerate", "[zonotope]") {
    Eigen::VectorXd c = Eigen::Vector2d::Zero();

    SECTION("full rank") {
        auto z = Zonotope::from_unit_box(2);
        REQUIRE_FALSE(z.degenerate());
    }

    SECTION("rank deficient") {
        Eigen::MatrixXd G(2, 1);
        G << 1, 0;
        Zonotope z(c, G);
        REQUIRE(z.degenerate());
    }
}

TEST_CASE("Zonotope bounded", "[zonotope]") {
    auto z = Zonotope::from_unit_box(3);
    REQUIRE(z.bounded());
}

TEST_CASE("Zonotope translate", "[zonotope]") {
    auto z = Zonotope::from_unit_box(2);
    Eigen::Vector2d v(1.0, 2.0);
    auto z2 = z.translate(v);
    REQUIRE(z2.center().isApprox(v));
    REQUIRE(z2.contains_point(v));
    REQUIRE_FALSE(z2.contains_point(Eigen::Vector2d(-1.0, -1.0)));
}

TEST_CASE("Zonotope linear transform", "[zonotope]") {
    auto z = Zonotope::from_unit_box(2);
    Eigen::MatrixXd M(2, 2);
    M << 2, 0, 0, 3;
    auto z2 = z.linear_transform(M);

    REQUIRE(z2.contains_point(Eigen::Vector2d(2.0, 3.0)));
    REQUIRE(z2.contains_point(Eigen::Vector2d(-2.0, -3.0)));
    REQUIRE_FALSE(z2.contains_point(Eigen::Vector2d(2.5, 0.0)));
}

TEST_CASE("Zonotope minkowski sum", "[zonotope]") {
    auto z1 = Zonotope::from_unit_box(2);
    auto z2 = Zonotope::from_unit_box(2);

    auto zsum = minkowski_sum(z1, z2);
    REQUIRE(zsum.n_generators() == 4);
    REQUIRE(zsum.contains_point(Eigen::Vector2d(2.0, 2.0)));
    REQUIRE(zsum.contains_point(Eigen::Vector2d(-2.0, -2.0)));
    REQUIRE_FALSE(zsum.contains_point(Eigen::Vector2d(2.5, 0.0)));
}

TEST_CASE("Zonotope intersection throws NotImplemented", "[zonotope]") {
    auto z1 = Zonotope::from_unit_box(2);
    auto z2 = Zonotope::from_unit_box(2);

    REQUIRE_THROWS_AS(intersection(z1, z2), NotImplemented);
    REQUIRE_THROWS_AS(do_intersect(z1, z2), NotImplemented);
}

TEST_CASE("Zonotope norm - unit box", "[zonotope][zonotope_norm]") {
    auto z = Zonotope::from_unit_box(2);

    SECTION("origin has norm 0") {
        REQUIRE(z.zonotope_norm(Eigen::Vector2d::Zero()) == Catch::Approx(0.0));
    }

    SECTION("boundary points have norm 1") {
        REQUIRE(z.zonotope_norm(Eigen::Vector2d(1.0, 0.0)) == Catch::Approx(1.0));
        REQUIRE(z.zonotope_norm(Eigen::Vector2d(0.0, 1.0)) == Catch::Approx(1.0));
        REQUIRE(z.zonotope_norm(Eigen::Vector2d(1.0, 1.0)) == Catch::Approx(1.0));
        REQUIRE(z.zonotope_norm(Eigen::Vector2d(-1.0, -1.0)) == Catch::Approx(1.0));
    }

    SECTION("interior points have norm < 1") {
        REQUIRE(z.zonotope_norm(Eigen::Vector2d(0.5, 0.3)) < 1.0);
    }

    SECTION("exterior points have norm > 1") {
        REQUIRE(z.zonotope_norm(Eigen::Vector2d(1.5, 0.0)) == Catch::Approx(1.5));
        REQUIRE(z.zonotope_norm(Eigen::Vector2d(2.0, 2.0)) == Catch::Approx(2.0));
    }
}

TEST_CASE("Zonotope norm - non-axis-aligned generators", "[zonotope][zonotope_norm]") {
    Eigen::VectorXd c = Eigen::Vector2d::Zero();
    Eigen::MatrixXd G(2, 2);
    G << 1, 1, 0, 1;
    Zonotope z(c, G);

    // Vertex [2,1] = G*[1,1], so norm should be 1
    REQUIRE(z.zonotope_norm(Eigen::Vector2d(2.0, 1.0)) == Catch::Approx(1.0));

    // [1, 0.5] = G*[0.5, 0.5], norm should be 0.5
    REQUIRE(z.zonotope_norm(Eigen::Vector2d(1.0, 0.5)) == Catch::Approx(0.5));

    // Unreachable direction (outside column space) returns infinity
    // For a full-rank G, all points are reachable, so test with rank-deficient G
    Eigen::MatrixXd G_degen(2, 1);
    G_degen << 1, 0;
    Zonotope z_degen(c, G_degen);
    REQUIRE(z_degen.zonotope_norm(Eigen::Vector2d(0.0, 1.0)) ==
            std::numeric_limits<double>::infinity());
}

TEST_CASE("Zonotope norm - no generators", "[zonotope][zonotope_norm]") {
    Eigen::VectorXd c = Eigen::Vector2d::Zero();
    Eigen::MatrixXd G(2, 0);
    Zonotope z(c, G);

    REQUIRE(z.zonotope_norm(Eigen::Vector2d::Zero()) == Catch::Approx(0.0));
    REQUIRE(z.zonotope_norm(Eigen::Vector2d(1.0, 0.0)) == std::numeric_limits<double>::infinity());
}

TEST_CASE("Zonotope to_vertices - unit box 2D", "[zonotope]") {
    auto z = Zonotope::from_unit_box(2);
    auto verts = z.to_vertices();

    // Unit box should have 4 vertices
    REQUIRE(verts.rows() == 4);
    REQUIRE(verts.cols() == 2);
}

TEST_CASE("Zonotope to_vertices - unit box 3D", "[zonotope]") {
    auto z = Zonotope::from_unit_box(3);
    auto verts = z.to_vertices();

    // 3D unit box should have 8 vertices
    REQUIRE(verts.rows() == 8);
    REQUIRE(verts.cols() == 3);
}
