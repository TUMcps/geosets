#include "catch2/catch_test_macros.hpp"
#include "geosets/binary_operations.h"
#include <Eigen/Dense>
#include <catch2/catch_all.hpp>
#include <geosets/exceptions.h>
#include <geosets/sets/interval.h>

using namespace geosets;

TEST_CASE("Interval returns both line intersection points", "[interval][line]") {
    const geosets::Interval interval = geosets::Interval::from_unit_box(2);

    const auto [lower, upper] =
        interval.line_intersection(Eigen::Vector2d::Zero(), Eigen::Vector2d(1.0, 0.0));

    REQUIRE(lower.isApprox(Eigen::Vector2d(-1.0, 0.0)));
    REQUIRE(upper.isApprox(Eigen::Vector2d(1.0, 0.0)));
}

TEST_CASE("Interval construction edge cases", "[interval][construction]") {
    SECTION("Custom interval construction") {
        Eigen::Vector2d lb(0, 1);
        Eigen::Vector2d ub(2, 3);

        auto interval = Interval(lb, ub);
        REQUIRE(interval.dim() == 2);
        REQUIRE(interval.lb.isApprox(lb));
        REQUIRE(interval.ub.isApprox(ub));
    }

    SECTION("Validation throws on invalid bounds") {
        Eigen::Vector2d lb(2, 1);
        Eigen::Vector2d ub(1, 3);

        REQUIRE_THROWS_AS(Interval(lb, ub), IntervalValidationException);
    }
}

TEST_CASE("Interval contains_point edge cases", "[interval][contains]") {
    auto interval = Interval::from_unit_box(2);

    SECTION("Boundary points are included") {
        Eigen::Vector2d point1(1.0, 0.0);
        Eigen::Vector2d point2(0.0, -1.0);
        Eigen::Vector2d point3(-1.0, 1.0);

        REQUIRE(interval.contains_point(point1));
        REQUIRE(interval.contains_point(point2));
        REQUIRE(interval.contains_point(point3));
    }

    SECTION("Points exactly on bounds are included") {
        Eigen::Vector2d lower_bound(-1.0, -1.0);
        Eigen::Vector2d upper_bound(1.0, 1.0);

        REQUIRE(interval.contains_point(lower_bound));
        REQUIRE(interval.contains_point(upper_bound));
    }
}

TEST_CASE("Interval support_function edge cases", "[interval][support]") {
    auto interval = Interval::from_unit_box(2);

    SECTION("Support vector non-uniqueness for axis-aligned directions") {
        auto support = interval.support_function(Eigen::Vector2d(1, 0));
        REQUIRE(support.vector.x() == 1.0);
        // y-coordinate is not unique - can be any value in [-1, 1]
        REQUIRE(support.vector.y() >= -1.0);
        REQUIRE(support.vector.y() <= 1.0);
        REQUIRE(support.value == Catch::Approx(1.0));
    }

    SECTION("Mixed direction support") {
        auto support = interval.support_function(Eigen::Vector2d(1, -1));
        REQUIRE(support.vector.isApprox(Eigen::Vector2d(1, -1)));
        REQUIRE(support.value == Catch::Approx(2.0)); // 1*1 + (-1)*(-1) = 2
    }
}

TEST_CASE("Interval vertex generation", "[interval][vertices]") {
    SECTION("2D interval generates correct number of vertices") {
        auto interval = Interval::from_unit_box(2);
        auto vertices = interval.to_vertices();

        REQUIRE(vertices.rows() == 4); // 2^2 = 4 vertices for 2D
        REQUIRE(vertices.cols() == 2);

        // Verify all vertices are within bounds
        for (int i = 0; i < vertices.rows(); ++i) {
            REQUIRE(vertices(i, 0) >= -1.0);
            REQUIRE(vertices(i, 0) <= 1.0);
            REQUIRE(vertices(i, 1) >= -1.0);
            REQUIRE(vertices(i, 1) <= 1.0);
        }
    }

    SECTION("3D interval generates correct number of vertices") {
        auto interval = Interval::from_unit_box(3);
        auto vertices = interval.to_vertices();

        REQUIRE(vertices.rows() == 8); // 2^3 = 8 vertices for 3D
        REQUIRE(vertices.cols() == 3);
    }
}

TEST_CASE("Interval matrix multiplication edge cases", "[interval][linear_transform]") {
    SECTION("Rotation preserves unit box structure") {
        auto interval = Interval::from_unit_box(2);
        // 90 degree rotation matrix
        Eigen::Matrix2d rotation;
        rotation << 0, -1, 1, 0;

        auto rotated = interval.linear_transform(rotation);
        // After 90 degree rotation, the interval should still be [-1,1]x[-1,1]
        REQUIRE(rotated.lb.isApprox(Eigen::Vector2d(-1, -1)));
        REQUIRE(rotated.ub.isApprox(Eigen::Vector2d(1, 1)));
    }
}

TEST_CASE("Interval binary operations", "[interval][binary_ops]") {
    SECTION("Intersection of overlapping intervals") {
        Interval a(Eigen::Vector2d(0, 0), Eigen::Vector2d(2, 2));
        Interval b(Eigen::Vector2d(1, 1), Eigen::Vector2d(3, 3));

        auto result = intersection(a, b);

        REQUIRE(result.lb.isApprox(Eigen::Vector2d(1, 1)));
        REQUIRE(result.ub.isApprox(Eigen::Vector2d(2, 2)));
    }

    SECTION("Intersection of identical intervals") {
        Interval a(Eigen::Vector2d(0, 0), Eigen::Vector2d(2, 2));
        Interval b(Eigen::Vector2d(0, 0), Eigen::Vector2d(2, 2));

        auto result = intersection(a, b);

        REQUIRE(result.lb.isApprox(Eigen::Vector2d(0, 0)));
        REQUIRE(result.ub.isApprox(Eigen::Vector2d(2, 2)));
    }

    SECTION("Intersection of non-overlapping intervals is empty") {
        Interval a(Eigen::Vector2d(0, 0), Eigen::Vector2d(1, 1));
        Interval b(Eigen::Vector2d(2, 2), Eigen::Vector2d(3, 3));

        REQUIRE(intersection(a, b).is_empty());
    }

    SECTION("Minkowski sum of intervals") {
        Interval a(Eigen::Vector2d(0, 0), Eigen::Vector2d(2, 2));
        Interval b(Eigen::Vector2d(1, 1), Eigen::Vector2d(3, 3));

        auto result = minkowski_sum(a, b);

        REQUIRE(result.lb.isApprox(Eigen::Vector2d(1, 1)));
        REQUIRE(result.ub.isApprox(Eigen::Vector2d(5, 5)));
    }

    SECTION("Minkowski sum with zero interval") {
        Interval a(Eigen::Vector2d(0, 0), Eigen::Vector2d(2, 2));
        Interval zero(Eigen::Vector2d(0, 0), Eigen::Vector2d(0, 0));

        auto result = minkowski_sum(a, zero);

        REQUIRE(result.lb.isApprox(Eigen::Vector2d(0, 0)));
        REQUIRE(result.ub.isApprox(Eigen::Vector2d(2, 2)));
    }

    SECTION("Minkowski sum is commutative") {
        Interval a(Eigen::Vector2d(1, 2), Eigen::Vector2d(3, 4));
        Interval b(Eigen::Vector2d(5, 6), Eigen::Vector2d(7, 8));

        auto sum_ab = minkowski_sum(a, b);
        auto sum_ba = minkowski_sum(b, a);

        REQUIRE(sum_ab.lb.isApprox(sum_ba.lb));
        REQUIRE(sum_ab.ub.isApprox(sum_ba.ub));
    }

    SECTION("Dimension mismatch throws") {
        Interval a2d(Eigen::Vector2d(0, 0), Eigen::Vector2d(1, 1));
        Interval b3d(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));

        REQUIRE_THROWS_AS(intersection(a2d, b3d), DimensionMismatchException);
        REQUIRE_THROWS_AS(minkowski_sum(a2d, b3d), DimensionMismatchException);
    }
}
