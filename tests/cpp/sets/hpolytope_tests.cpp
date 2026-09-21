#include "catch2/catch_test_macros.hpp"
#include <Eigen/Dense>
#include <catch2/catch_all.hpp>
#include <geosets/linalg_utils.h>
#include <geosets/sets/hpolytope.h>

using namespace geosets;

TEST_CASE("HPolytope returns both line intersection points", "[hpolytope][line]") {
    const HPolytope box = HPolytope::from_unit_box(2);

    const auto [lower, upper] =
        box.line_intersection(Eigen::Vector2d::Zero(), Eigen::Vector2d(1.0, 0.0));

    REQUIRE(lower.isApprox(Eigen::Vector2d(-1.0, 0.0)));
    REQUIRE(upper.isApprox(Eigen::Vector2d(1.0, 0.0)));
}

TEST_CASE("Fourier-Motzkin elimination dimension reduction", "[hpolytope][fourier_motzkin]") {
    auto box = HPolytope::from_unit_box(3);

    auto [A_proj, b_proj] = linalg::fourier_motzkin_elimination(box.A, box.b, 1);

    REQUIRE(A_proj.cols() == 2);
    REQUIRE(A_proj.rows() == b_proj.size());

    HPolytope projected(A_proj, b_proj);
    REQUIRE_FALSE(projected.is_empty());
}

TEST_CASE("Fourier-Motzkin with negative bounds", "[hpolytope][fourier_motzkin]") {
    // Box: -2 <= x <= 0, -3 <= y <= -1
    Eigen::MatrixXd A(4, 2);
    A << 1, 0, // x <= 0
        0, 1,  // y <= -1
        -1, 0, // x >= -2
        0, -1; // y >= -3
    Eigen::VectorXd b(4);
    b << 0, -1, 2, 3;

    auto [A_proj, b_proj] = linalg::fourier_motzkin_elimination(A, b, 0);

    REQUIRE(A_proj.cols() == 1);
    HPolytope projected(A_proj, b_proj);
    REQUIRE_FALSE(projected.is_empty());
}

TEST_CASE("Fourier-Motzkin with zero coefficients", "[hpolytope][fourier_motzkin]") {
    // Some constraints independent of elimination dimension
    Eigen::MatrixXd A(3, 2);
    A << 1, 0, 0, 1, 0, -1;
    Eigen::VectorXd b = Eigen::VectorXd::Ones(3);

    auto [A_proj, b_proj] = linalg::fourier_motzkin_elimination(A, b, 0);

    REQUIRE(A_proj.cols() == 1);
    REQUIRE(A_proj.rows() >= 2);
}

TEST_CASE("Fourier-Motzkin all positive/negative coefficients", "[hpolytope][fourier_motzkin]") {
    SECTION("All positive") {
        Eigen::MatrixXd A(2, 2);
        A << 1, 1, 1, 0;
        Eigen::VectorXd b = Eigen::Vector2d(2, 1);

        auto [A_proj, b_proj] = linalg::fourier_motzkin_elimination(A, b, 0);
        REQUIRE(A_proj.cols() == 1);
    }

    SECTION("All negative") {
        Eigen::MatrixXd A(2, 2);
        A << -1, 1, -1, 0;
        Eigen::VectorXd b = Eigen::Vector2d(2, 1);

        auto [A_proj, b_proj] = linalg::fourier_motzkin_elimination(A, b, 0);
        REQUIRE(A_proj.cols() == 1);
    }
}

TEST_CASE("Linear transform dimension reduction 3D to 2D", "[hpolytope][linear_transform]") {
    auto box = HPolytope::from_unit_box(3);

    // Project onto first two dimensions
    Eigen::MatrixXd projection = Eigen::MatrixXd::Zero(2, 3);
    projection(0, 0) = 1.0;
    projection(1, 1) = 1.0;

    box.linear_transform_(projection);

    REQUIRE(box.dim() == 2);
    REQUIRE_FALSE(box.is_empty());
    REQUIRE(box.bounded());

    auto center = box.center();
    REQUIRE(center.isApprox(Eigen::Vector2d(0.0, 0.0)));
}

TEST_CASE("Linear transform dimension reduction 3D to 1D", "[hpolytope][linear_transform]") {
    auto box = HPolytope::from_unit_box(3);

    // Project onto first dimension only
    Eigen::MatrixXd projection = Eigen::MatrixXd::Zero(1, 3);
    projection(0, 0) = 1.0;

    box.linear_transform_(projection);

    REQUIRE(box.dim() == 1);
    REQUIRE_FALSE(box.is_empty());
    REQUIRE(box.bounded());
}

TEST_CASE("Linear transform with negative bounds", "[hpolytope][linear_transform]") {
    // Create a box [-2, -1] x [-3, -2] x [0, 1] in 3D
    Eigen::MatrixXd A(6, 3);
    A << 1, 0, 0, // x <= -1
        0, 1, 0,  // y <= -2
        0, 0, 1,  // z <= 1
        -1, 0, 0, // x >= -2
        0, -1, 0, // y >= -3
        0, 0, -1; // z >= 0
    Eigen::VectorXd b(6);
    b << -1, -2, 1, 2, 3, 0;

    HPolytope poly(A, b);

    // Project to first two dimensions
    Eigen::MatrixXd projection = Eigen::MatrixXd::Zero(2, 3);
    projection(0, 0) = 1.0;
    projection(1, 1) = 1.0;

    poly.linear_transform_(projection);

    REQUIRE(poly.dim() == 2);
    REQUIRE_FALSE(poly.is_empty());

    // Center should be at [-1.5, -2.5]
    auto center = poly.center();
    REQUIRE(center.isApprox(Eigen::Vector2d(-1.5, -2.5), 1e-6));
}

TEST_CASE("Linear transform with scaling and projection", "[hpolytope][linear_transform]") {
    auto box = HPolytope::from_unit_box(3);

    // Project onto first two dimensions with scaling
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(2, 3);
    M(0, 0) = 2.0; // Scale x by 2
    M(1, 1) = 0.5; // Scale y by 0.5

    box.linear_transform_(M);

    REQUIRE(box.dim() == 2);
    REQUIRE_FALSE(box.is_empty());
    REQUIRE(box.bounded());

    // Check support function values to verify the projected box has correct bounds
    auto support_x = box.support_function(Eigen::Vector2d(1.0, 0.0));
    auto support_y = box.support_function(Eigen::Vector2d(0.0, 1.0));

    REQUIRE(std::abs(support_x.value - 2.0) < 1e-6); // Max x should be 2.0
    REQUIRE(std::abs(support_y.value - 0.5) < 1e-6); // Max y should be 0.5
}

TEST_CASE("Linear transform with rotation and projection", "[hpolytope][linear_transform]") {
    auto box = HPolytope::from_unit_box(3);

    // Combine rotation in xy-plane with projection
    const double angle = M_PI / 4.0;
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(2, 3);
    M(0, 0) = std::cos(angle);
    M(0, 1) = -std::sin(angle);
    M(1, 0) = std::sin(angle);
    M(1, 1) = std::cos(angle);

    box.linear_transform_(M);

    REQUIRE(box.dim() == 2);
    REQUIRE_FALSE(box.is_empty());
    REQUIRE(box.bounded());
}
