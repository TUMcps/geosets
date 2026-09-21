#include "../utils.h"
#include "catch2/catch_test_macros.hpp"
#include "geosets/binary_operations.h"
#include <Eigen/Dense>
#include <catch2/catch_all.hpp>
#include <geosets/sets/polygon.h>

using namespace geosets;

constexpr double tol = 1e-9;

TEST_CASE("Polygon unary operations", "[polygon][unary_ops]") {
    SECTION("Volume and degenerate checks") {
        auto box = Polygon::from_unit_box(2);
        REQUIRE(std::abs(box.volume() - 4.0) < tol);
        REQUIRE_FALSE(box.degenerate());

        Eigen::MatrixXd line_segment(4, 2);
        line_segment << -1.0, -1.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0;
        Polygon degenerate(line_segment, false);
        REQUIRE(degenerate.degenerate());
        REQUIRE(degenerate.volume() == 0.0);
    }

    SECTION("Containment supports tolerance near boundary") {
        auto box = Polygon::from_unit_box(2);

        REQUIRE(box.contains_point(Eigen::Vector2d(-0.5, 0.25)));
        REQUIRE(box.contains_point(Eigen::Vector2d(-1.0, -1.0)));
        REQUIRE_FALSE(box.contains_point(Eigen::Vector2d(-1.0 - 1e-6, 0.0)));
        REQUIRE(box.contains_point(Eigen::Vector2d(-1.0 - 1e-10, 0.0), 1e-9));
    }

    SECTION("Support function with negative coordinates") {
        Eigen::MatrixXd vertices(5, 2);
        vertices << -2.0, -1.0, 1.0, -1.0, 1.0, 2.0, -2.0, 2.0, -2.0, -1.0;
        Polygon polygon(vertices);

        auto support = polygon.support_function(Eigen::Vector2d(-1.0, -1.0));
        REQUIRE(support.vector.isApprox(Eigen::Vector2d(-2.0, -1.0)));
        REQUIRE(std::abs(support.value - 3.0) < tol);
    }

    SECTION("Linear transform with reflection") {
        auto box = Polygon::from_unit_box(2);
        Eigen::Matrix2d reflection;
        reflection << -1.0, 0.0, 0.0, 1.0;

        auto transformed = box.linear_transform(reflection);
        REQUIRE(transformed.center().isApprox(Eigen::Vector2d::Zero()));
        REQUIRE(std::abs(transformed.volume() - 4.0) < tol);
    }
}

TEST_CASE("Polygon binary operations", "[polygon][binary_ops]") {
    SECTION("Intersection of translated boxes") {
        auto a = Polygon::from_unit_box(2).translate(Eigen::Vector2d(0.5, 0.5));
        auto b = Polygon::from_unit_box(2).translate(Eigen::Vector2d(-0.5, -0.5));
        auto result = intersection(a, b);

        Eigen::MatrixXd expected(4, 2);
        expected << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;

        REQUIRE(vertices_to_set(result.to_vertices()) == vertices_to_set(expected));
        REQUIRE(result.center().isApprox(Eigen::Vector2d::Zero()));
        REQUIRE(std::abs(result.volume() - 1.0) < tol);
    }

    SECTION("Intersection of disjoint polygons is empty") {
        auto a = Polygon::from_unit_box(2).translate(Eigen::Vector2d(-3.0, 0.0));
        auto b = Polygon::from_unit_box(2).translate(Eigen::Vector2d(3.0, 0.0));
        REQUIRE(intersection(a, b).is_empty());
    }

    SECTION("Minkowski sum of boxes") {
        auto a = Polygon::from_unit_box(2);
        auto b = Polygon::from_unit_box(2);
        auto sum = minkowski_sum(a, b);

        Eigen::MatrixXd expected(4, 2);
        expected << 2.0, -2.0, 2.0, 2.0, -2.0, 2.0, -2.0, -2.0;

        REQUIRE(vertices_to_set(sum.to_vertices()) == vertices_to_set(expected));
        REQUIRE(sum.center().isApprox(Eigen::Vector2d::Zero()));
        REQUIRE(std::abs(sum.volume() - 16.0) < tol);
    }

    SECTION("Set difference of overlapping boxes") {
        auto a = Polygon::from_unit_box(2);
        auto b = Polygon::from_unit_box(2).translate(Eigen::Vector2d(1.0, 0.0));
        auto result = set_difference(a, b);

        Eigen::MatrixXd expected(4, 2);
        expected << -1.0, -1.0, 0.0, -1.0, 0.0, 1.0, -1.0, 1.0;

        REQUIRE(vertices_to_set(result.to_vertices()) == vertices_to_set(expected));
        REQUIRE(result.center().isApprox(Eigen::Vector2d(-0.5, 0.0)));
        REQUIRE(std::abs(result.volume() - 2.0) < tol);
    }

    SECTION("Set difference of disjoint polygons returns minuend") {
        auto a = Polygon::from_unit_box(2).translate(Eigen::Vector2d(-2.0, -1.0));
        auto b = Polygon::from_unit_box(2).translate(Eigen::Vector2d(3.0, 2.0));
        auto result = set_difference(a, b);

        REQUIRE(vertices_to_set(result.to_vertices()) == vertices_to_set(a.to_vertices()));
        REQUIRE(result.center().isApprox(a.center()));
        REQUIRE(std::abs(result.volume() - a.volume()) < tol);
    }

    SECTION("Set difference that becomes disconnected returns first component") {
        auto a = Polygon::from_unit_box(2);
        Eigen::MatrixXd strip_vertices(4, 2);
        strip_vertices << -0.25, -2.0, 0.25, -2.0, 0.25, 2.0, -0.25, 2.0;
        Polygon strip(strip_vertices);

        auto result = set_difference(a, strip);
        REQUIRE_FALSE(result.is_empty());
        REQUIRE(result.volume() > 0.0);
        REQUIRE(result.volume() < a.volume());
    }

    SECTION("Set difference with complete overlap returns empty polygon") {
        auto a = Polygon::from_unit_box(2);
        auto result = set_difference(a, a);

        REQUIRE(result.is_empty());
    }
}
