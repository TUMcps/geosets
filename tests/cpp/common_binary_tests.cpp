#include "catch2/catch_template_test_macros.hpp"
#include "catch2/catch_test_macros.hpp"
#include "definitions.h"
#include "geosets/sets/hpolytope.h"
#include "geosets/sets/polygon.h"
#include "utils.h"
#include <Eigen/Dense>
#include <catch2/catch_all.hpp>
#include <geosets/binary_operations.h>
#include <geosets/sets/interval.h>
#include <geosets/sets/vpolytope.h>
#include <geosets/sets/zonotope.h>

// Test types that support binary operations
using BinaryOpSetTypes = std::tuple<geosets::Interval, geosets::VPolytope, geosets::HPolytope,
                                    geosets::Polygon, geosets::Zonotope>;

TEMPLATE_LIST_TEST_CASE("Binary operations - Dimension compatibility", "[binary_ops][dimensions]",
                        BinaryOpSetTypes) {
    if (std::is_same<TestType, geosets::Polygon>::value) {
        return;
    }

    TestType a2d = TestType::from_unit_box(2);
    TestType b3d = TestType::from_unit_box(3);

    REQUIRE_THROWS_AS_OR_SKIP(intersection(a2d, b3d), geosets::DimensionMismatchException);
    REQUIRE_THROWS_AS_OR_SKIP(do_intersect(a2d, b3d), geosets::DimensionMismatchException);
    REQUIRE_THROWS_AS_OR_SKIP(minkowski_sum(a2d, b3d), geosets::DimensionMismatchException);
}

TEMPLATE_LIST_TEST_CASE("Binary operations - Intersect predicate", "[binary_ops][intersect]",
                        BinaryOpSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        TestType a = TestType::from_unit_box(dim);
        TestType b = TestType::from_unit_box(dim).translate(Eigen::VectorXd::Constant(dim, 0.5));
        TestType disjoint =
            TestType::from_unit_box(dim).translate(Eigen::VectorXd::Constant(dim, 3.0));

        Eigen::VectorXd touching_translation = Eigen::VectorXd::Zero(dim);
        touching_translation(0) = 2.0;
        TestType touching = TestType::from_unit_box(dim).translate(touching_translation);

        REQUIRE_OR_SKIP(do_intersect(a, b));
        REQUIRE_OR_SKIP(do_intersect(b, a));
        REQUIRE_FALSE_OR_SKIP(do_intersect(a, disjoint));
        REQUIRE_FALSE_OR_SKIP(do_intersect(disjoint, a));
        REQUIRE_OR_SKIP(do_intersect(a, touching));
    }
}

TEMPLATE_LIST_TEST_CASE("Binary operations - Intersection", "[binary_ops][intersection]",
                        BinaryOpSetTypes) {
    SECTION("Intersection of identical sets returns the same set") {
        auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
        for (size_t dim = dim_start; dim < dim_end; ++dim) {
            auto a = TestType::from_unit_box(dim);
            skip_if_unimplemented([&] {
                auto result = intersection(a, a);
                REQUIRE(result.dim() == a.dim());
                REQUIRE(result.center() == a.center());
                REQUIRE(std::abs(result.volume() - a.volume()) <= tol_);
                REQUIRE_FALSE(result.is_empty());
            });
        }
    }

    SECTION("Intersection result is contained in both operands") {
        auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
        for (size_t dim = dim_start; dim < dim_end; ++dim) {
            TestType a = TestType::from_unit_box(dim);
            TestType b =
                TestType::from_unit_box(dim).translate(Eigen::VectorXd(dim).setConstant(0.5));

            skip_if_unimplemented([&] {
                auto result = intersection(a, b);
                REQUIRE(a.contains_point(result.center()));
                REQUIRE(b.contains_point(result.center()));

                // Should be commutative
                auto ab = intersection(a, b);
                auto ba = intersection(b, a);
                REQUIRE(ab.center() == ba.center());
                REQUIRE(ab.volume() == ba.volume());
            });
        }
    }
}

TEMPLATE_LIST_TEST_CASE("Binary operations - Minkowski sum", "[binary_ops][minkowski_sum]",
                        BinaryOpSetTypes) {

    SECTION("2D unit boxes") {
        skip_if_unimplemented([&] {
            auto a = TestType::from_unit_box(2);
            auto b = TestType::from_unit_box(2);

            auto ab = minkowski_sum(a, b);
            auto ba = minkowski_sum(b, a);

            Eigen::MatrixXd expected(4, 2);
            expected << 2.0, -2.0, 2.0, 2.0, -2.0, 2.0, -2.0, -2.0;

            REQUIRE(vertices_to_set(ab.to_vertices()) == vertices_to_set(expected));
            REQUIRE(vertices_to_set(ba.to_vertices()) == vertices_to_set(expected));

            REQUIRE(ab.center() == Eigen::Vector2d::Zero(2));
            REQUIRE(ba.center() == Eigen::Vector2d::Zero(2));
            REQUIRE(std::abs(ab.volume() - 16.0) < tol_);
            REQUIRE(std::abs(ba.volume() - 16.0) < tol_);
        });
    }

    if (std::is_same<TestType, geosets::Polygon>::value) {
        return;
    }

    SECTION("3D unit boxes") {
        skip_if_unimplemented([&] {
            auto a = TestType::from_unit_box(3);
            auto b = TestType::from_unit_box(3);

            auto ab = minkowski_sum(a, b);
            auto ba = minkowski_sum(b, a);

            Eigen::MatrixXd expected(8, 3);
            expected << 2.0, -2.0, -2.0, 2.0, -2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, -2.0, -2.0, 2.0,
                2.0, -2.0, 2.0, -2.0, -2.0, -2.0, 2.0, -2.0, -2.0, -2.0;

            REQUIRE(vertices_to_set(ab.to_vertices()) == vertices_to_set(expected));
            REQUIRE(vertices_to_set(ba.to_vertices()) == vertices_to_set(expected));
            REQUIRE(ab.center() == Eigen::Vector3d::Zero(3));
            REQUIRE(ba.center() == Eigen::Vector3d::Zero(3));
            REQUIRE(std::abs(ab.volume() - 64.0) < tol_);
            REQUIRE(std::abs(ba.volume() - 64.0) < tol_);
        });
    }
}

TEST_CASE("Minkowski Sum with random feasible", "[binary_ops][minkowski_sum]") {
    SECTION("HPolytope") {
        for (size_t d = 2, n = 3; d <= 3; ++d, n *= 2) {
            auto a = geosets::HPolytope::from_random(d, n);
            auto b = geosets::HPolytope::from_random(d, n);
            auto sum = minkowski_sum(a, b);

            REQUIRE_FALSE(sum.is_empty());
            REQUIRE_FALSE(sum.degenerate());
            REQUIRE(sum.volume() > 0.0);
        }
    }
}

TEMPLATE_LIST_TEST_CASE("Binary operations - Set difference", "[binary_ops][set_difference]",
                        BinaryOpSetTypes) {
    if constexpr (std::is_same_v<TestType, geosets::Polygon>) {
        SECTION("Overlapping boxes") {
            auto a = TestType::from_unit_box(2);
            auto b = TestType::from_unit_box(2).translate(Eigen::Vector2d(1.0, 0.0));
            auto result = set_difference(a, b);

            Eigen::MatrixXd expected(4, 2);
            expected << -1.0, -1.0, 0.0, -1.0, 0.0, 1.0, -1.0, 1.0;

            REQUIRE(vertices_to_set(result.to_vertices()) == vertices_to_set(expected));
            REQUIRE(result.center().isApprox(Eigen::Vector2d(-0.5, 0.0)));
            REQUIRE(std::abs(result.volume() - 2.0) < tol_);
        }

        SECTION("Disjoint sets return minuend") {
            auto a = TestType::from_unit_box(2).translate(Eigen::Vector2d(-2.0, -1.0));
            auto b = TestType::from_unit_box(2).translate(Eigen::Vector2d(3.0, 2.0));
            auto result = set_difference(a, b);

            REQUIRE(vertices_to_set(result.to_vertices()) == vertices_to_set(a.to_vertices()));
            REQUIRE(result.center().isApprox(a.center()));
            REQUIRE(std::abs(result.volume() - a.volume()) < tol_);
        }

        SECTION("Complete overlap returns empty") {
            auto a = TestType::from_unit_box(2);
            auto result = set_difference(a, a);
            REQUIRE(result.is_empty());
        }

        SECTION("Disconnected result returns first component") {
            auto a = TestType::from_unit_box(2);
            Eigen::MatrixXd strip_vertices(4, 2);
            strip_vertices << -0.25, -2.0, 0.25, -2.0, 0.25, 2.0, -0.25, 2.0;
            TestType strip(strip_vertices);
            auto result = set_difference(a, strip);

            REQUIRE_FALSE(result.is_empty());
            REQUIRE(result.volume() > 0.0);
            REQUIRE(result.volume() < a.volume());
        }
    }
}
