#include "catch2/catch_template_test_macros.hpp"
#include "catch2/catch_test_macros.hpp"
#include "definitions.h"
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <catch2/catch_all.hpp>
#include <geosets/sets/hpolytope.h>
#include <geosets/sets/interval.h>
#include <geosets/sets/polygon.h>
#include <geosets/sets/vpolytope.h>

TEMPLATE_LIST_TEST_CASE("Geoset dimension matches unit box construction", "[geosets][dim]",
                        GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        auto set = TestType::from_unit_box(dim);
        REQUIRE(set.dim() == dim);
    }
}

TEMPLATE_LIST_TEST_CASE("Empty", "[unary_ops][empty]", GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        auto set = TestType::from_unit_box(dim);
        REQUIRE_FALSE(set.is_empty());
    }
}

TEMPLATE_LIST_TEST_CASE("Empty instantiation", "[unary_ops][empty][instantiation]", GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        auto empty_set = TestType::make_empty(dim);
        REQUIRE(empty_set.dim() == dim);
        REQUIRE(empty_set.is_empty());
        REQUIRE_FALSE(empty_set.is_infinite());
    }
}

TEMPLATE_LIST_TEST_CASE("Empty set properties", "[unary_ops][empty][properties]", GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        auto empty_set = TestType::make_empty(dim);

        // Empty sets are bounded
        REQUIRE(empty_set.bounded());

        // Empty sets have zero volume
        skip_if_unimplemented([&] { REQUIRE(empty_set.volume() == 0.0); });

        // Empty sets don't contain any points
        skip_if_unimplemented([&] {
            auto test_point = Eigen::VectorXd::Zero(dim);
            REQUIRE_FALSE(empty_set.contains_point(test_point));
        });
    }
}

TEMPLATE_LIST_TEST_CASE("Empty set operations", "[unary_ops][empty][operations]", GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        auto empty_set = TestType::make_empty(dim);

        // Translation doesn't change emptiness
        auto translated = empty_set.translate(Eigen::VectorXd::Ones(dim));
        REQUIRE(translated.is_empty());

        // Linear transform doesn't change emptiness
        skip_if_unimplemented([&] {
            auto transformed =
                empty_set.linear_transform(Eigen::MatrixXd::Identity(dim, dim) * 2.0);
            REQUIRE(transformed.is_empty());
        });

        REQUIRE(empty_set.support_function(Eigen::VectorXd::Ones(dim)).value ==
                -std::numeric_limits<double>::infinity());
        REQUIRE_THROWS(empty_set.center());
        REQUIRE_THROWS(empty_set.to_vertices());
    }
}

TEMPLATE_LIST_TEST_CASE("Bounded", "[unary_ops][bounded]", GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        auto set = TestType::from_unit_box(dim);
        REQUIRE(set.bounded());
    }
}

TEMPLATE_LIST_TEST_CASE("Infinite instantiation", "[unary_ops][infinite][instantiation]",
                        GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        auto infinite_set = TestType::make_infinite(dim);
        REQUIRE(infinite_set.dim() == dim);
        REQUIRE(infinite_set.is_infinite());
        REQUIRE_FALSE(infinite_set.is_empty());
    }
}

TEMPLATE_LIST_TEST_CASE("Infinite set properties", "[unary_ops][infinite][properties]",
                        GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        auto infinite_set = TestType::make_infinite(dim);

        // Infinite sets are unbounded
        REQUIRE_FALSE(infinite_set.bounded());

        // Infinite sets have infinite volume
        skip_if_unimplemented(
            [&] { REQUIRE(infinite_set.volume() == std::numeric_limits<double>::infinity()); });

        // Infinite sets contain all points
        skip_if_unimplemented([&] {
            for (int i = 0; i < 5; ++i) {
                auto test_point = Eigen::VectorXd::Random(dim) * 1000.0; // Large random points
                REQUIRE(infinite_set.contains_point(test_point));
            }
        });

        // Infinite sets are not degenerate
        skip_if_unimplemented([&] { REQUIRE_FALSE(infinite_set.degenerate()); });
    }
}

TEMPLATE_LIST_TEST_CASE("Infinite set operations", "[unary_ops][infinite][operations]",
                        GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        auto infinite_set = TestType::make_infinite(dim);

        // Translation doesn't change infiniteness
        auto translated = infinite_set.translate(Eigen::VectorXd::Ones(dim));
        REQUIRE(translated.is_infinite());

        // Linear transform doesn't change infiniteness
        skip_if_unimplemented([&] {
            auto transformed =
                infinite_set.linear_transform(Eigen::MatrixXd::Identity(dim, dim) * 2.0);
            REQUIRE(transformed.is_infinite());
        });

        // Center, support function, and to_vertices should throw for infinite sets
        REQUIRE_THROWS(infinite_set.center());
        REQUIRE_THROWS(infinite_set.to_vertices());

        // Support function returns infinite value
        auto support = infinite_set.support_function(Eigen::VectorXd::Ones(dim));
        REQUIRE(support.value == std::numeric_limits<double>::infinity());
    }
}

TEMPLATE_LIST_TEST_CASE("Center", "[unary_ops][center]", GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        auto set = TestType::from_unit_box(dim);
        REQUIRE(set.center().isApprox(Eigen::VectorXd::Zero(dim)));
    }
}

TEMPLATE_LIST_TEST_CASE("Volume", "[unary_ops][volume]", GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        auto set = TestType::from_unit_box(dim);
        skip_if_unimplemented([&] {
            auto volume = set.volume();
            REQUIRE(std::abs(volume - std::pow(2.0, dim)) < tol_);
        });
    }
}

TEMPLATE_LIST_TEST_CASE("Translate and center", "[unary_ops][translate]", GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        auto set = TestType::from_unit_box(dim);
        REQUIRE(set.center().isApprox(Eigen::VectorXd::Zero(dim)));

        for (int i = 0; i < 10; ++i) {
            auto vector = Eigen::VectorXd::Random(dim).eval();
            auto translated = set.translate(vector);
            auto center = translated.center();

            REQUIRE(translated.center().isApprox(vector));
        }
    }
}

TEMPLATE_LIST_TEST_CASE("Degenerate", "[unary_ops][degenerate]", GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        auto set = TestType::from_unit_box(dim);
        skip_if_unimplemented([&] {
            // Needs to be split in order for the skip_if_unimplemented to work!
            bool is_degenerate = set.degenerate();
            REQUIRE_FALSE(is_degenerate);
        });
    }
}

TEMPLATE_LIST_TEST_CASE("Linear transform", "[unary_ops][linear_transform]", GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        auto set = TestType::from_unit_box(dim);
        skip_if_unimplemented([&] {
            set.linear_transform_(Eigen::MatrixXd::Identity(dim, dim) * 2.0);
            auto center = set.center();
            REQUIRE(center.isApprox(Eigen::VectorXd::Zero(dim)));
        });
    }
}

TEMPLATE_LIST_TEST_CASE("Linear transform - dimension reduction",
                        "[unary_ops][linear_transform][dimensions]", GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        if (dim < 3u) continue; // Need at least 3D for projection to 2D

        auto set = TestType::from_unit_box(dim);
        skip_if_unimplemented([&] {
            // Project to first two dimensions
            Eigen::MatrixXd projection = Eigen::MatrixXd::Zero(2, dim);
            projection(0, 0) = 1.0;
            projection(1, 1) = 1.0;

            set.linear_transform_(projection);
            REQUIRE(set.dim() == 2);
            REQUIRE_FALSE(set.is_empty());

            // Center should still be at origin
            auto center = set.center();
            REQUIRE(center.isApprox(Eigen::VectorXd::Zero(2)));
        });
    }
}

TEMPLATE_LIST_TEST_CASE("Support function", "[unary_ops][support_function]", GeoSetTypes) {
    auto set = TestType::from_unit_box(2);

    // [1, 0] vector
    skip_if_unimplemented([&] {
        {
            auto support = set.support_function(Eigen::Vector2d(1.0, 0.0));
            REQUIRE(support.vector.x() == 1.0);
            // y-coordinate is not unique
            REQUIRE(support.vector.y() >= -1.0);
            REQUIRE(support.vector.y() <= 1.0);
            REQUIRE(support.value == 1.0);
        }

        {
            auto support = set.support_function(Eigen::Vector2d(0.0, 1.0));
            REQUIRE(support.vector.y() == 1.0);
            // x-coordinate is not unique
            REQUIRE(support.vector.x() >= -1.0);
            REQUIRE(support.vector.x() <= 1.0);
            REQUIRE(support.value == 1.0);
        }

        {
            auto support = set.support_function(Eigen::Vector2d(1.0, 1.0));
            CHECK(support.vector.isApprox(Eigen::Vector2d(1.0, 1.0)));
            CHECK(support.value == 2.0);
        }
    });
}

TEMPLATE_LIST_TEST_CASE("Boundary point", "[unary_ops][boundary_point]", GeoSetTypes) {
    auto set = TestType::from_unit_box(2);

    skip_if_unimplemented([&] {
        // Axis-aligned direction: boundary of [-1,1]^2 along x-axis from origin
        auto bp = set.boundary_point(Eigen::Vector2d(1.0, 0.0));
        REQUIRE(bp.isApprox(Eigen::Vector2d(1.0, 0.0), 1e-6));

        // Negative direction
        bp = set.boundary_point(Eigen::Vector2d(-1.0, 0.0));
        REQUIRE(bp.isApprox(Eigen::Vector2d(-1.0, 0.0), 1e-6));

        // Diagonal: ray from origin along (1,1) hits boundary of [-1,1]^2 at (1,1)
        bp = set.boundary_point(Eigen::Vector2d(1.0, 1.0));
        REQUIRE(bp.isApprox(Eigen::Vector2d(1.0, 1.0), 1e-6));

        // Non-unit direction: (2,0) should still give boundary at (1,0)
        bp = set.boundary_point(Eigen::Vector2d(2.0, 0.0));
        REQUIRE(bp.isApprox(Eigen::Vector2d(1.0, 0.0), 1e-6));

        // Boundary point must be contained in the set
        REQUIRE(set.contains_point(bp));
    });
}

TEMPLATE_LIST_TEST_CASE("Boundary point - random directions on unit box",
                        "[unary_ops][boundary_point]", GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        auto set = TestType::from_unit_box(dim);

        skip_if_unimplemented([&] {
            for (int i = 0; i < 20; ++i) {
                Eigen::VectorXd dir = Eigen::VectorXd::Random(dim);
                if (dir.norm() < 1e-9) continue;

                auto bp = set.boundary_point(dir);
                REQUIRE(bp.template lpNorm<Eigen::Infinity>() == Catch::Approx(1.0).margin(1e-6));
                REQUIRE(set.contains_point(bp));
            }
        });
    }
}

TEMPLATE_LIST_TEST_CASE("Boundary point - translated set", "[unary_ops][boundary_point]",
                        GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        auto set = TestType::from_unit_box(dim);
        Eigen::VectorXd offset = Eigen::VectorXd::Constant(dim, 2.0);
        auto translated = set.translate(offset);

        skip_if_unimplemented([&] {
            // Boundary along first axis from translated center
            Eigen::VectorXd dir = Eigen::VectorXd::Zero(dim);
            dir(0) = 1.0;

            auto bp = translated.boundary_point(dir);

            // Should be at center + (1,0,...,0) = (3, 2, ..., 2)
            Eigen::VectorXd expected = offset;
            expected(0) += 1.0;
            REQUIRE(bp.isApprox(expected, 1e-6));
            REQUIRE(translated.contains_point(bp));
        });
    }
}

TEMPLATE_LIST_TEST_CASE("Boundary point - empty and infinite", "[unary_ops][boundary_point]",
                        GeoSetTypes) {
    auto empty = TestType::make_empty(2);
    REQUIRE_THROWS(empty.boundary_point(Eigen::Vector2d(1.0, 0.0)));

    auto inf = TestType::make_infinite(2);
    REQUIRE_THROWS(inf.boundary_point(Eigen::Vector2d(1.0, 0.0)));
}

TEMPLATE_LIST_TEST_CASE("Tangent hyperplane - axis-aligned directions",
                        "[unary_ops][tangent_hyperplane]", GeoSetTypes) {
    auto set = TestType::from_unit_box(2);

    skip_if_unimplemented([&] {
        // Direction (1, 0): should hit upper bound on x, plane x = 1
        {
            auto hp = set.tangent_hyperplane(Eigen::Vector2d(1.0, 0.0));
            REQUIRE(hp.a.isApprox(Eigen::Vector2d(1.0, 0.0)));
            REQUIRE(hp.b == Catch::Approx(1.0).margin(tol_));
        }

        // Direction (-1, 0): should hit lower bound on x, plane -x = 1
        {
            auto hp = set.tangent_hyperplane(Eigen::Vector2d(-1.0, 0.0));
            REQUIRE(hp.a.isApprox(Eigen::Vector2d(-1.0, 0.0)));
            REQUIRE(hp.b == Catch::Approx(1.0).margin(tol_));
        }

        // Direction (0, 1): should hit upper bound on y, plane y = 1
        {
            auto hp = set.tangent_hyperplane(Eigen::Vector2d(0.0, 1.0));
            REQUIRE(hp.a.isApprox(Eigen::Vector2d(0.0, 1.0)));
            REQUIRE(hp.b == Catch::Approx(1.0).margin(tol_));
        }
    });
}

TEMPLATE_LIST_TEST_CASE("Tangent hyperplane - support point on plane",
                        "[unary_ops][tangent_hyperplane]", GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        auto set = TestType::from_unit_box(dim);

        skip_if_unimplemented([&] {
            for (int i = 0; i < 20; ++i) {
                Eigen::VectorXd dir = Eigen::VectorXd::Random(dim);
                if (dir.norm() < 1e-9) continue;

                auto hp = set.tangent_hyperplane(dir);
                auto sp = set.support_function(dir);

                // Support point must lie on the tangent hyperplane: a^T x = b
                REQUIRE(hp.a.dot(sp.vector) == Catch::Approx(hp.b).margin(1e-6));

                // Normal must point outward (same halfspace as direction)
                REQUIRE(hp.a.dot(dir) > -1e-9);
            }
        });
    }
}

TEMPLATE_LIST_TEST_CASE("Tangent hyperplane - empty and infinite",
                        "[unary_ops][tangent_hyperplane]", GeoSetTypes) {
    auto empty = TestType::make_empty(2);
    REQUIRE_THROWS(empty.tangent_hyperplane(Eigen::Vector2d(1.0, 0.0)));

    auto inf = TestType::make_infinite(2);
    REQUIRE_THROWS(inf.tangent_hyperplane(Eigen::Vector2d(1.0, 0.0)));
}

TEMPLATE_LIST_TEST_CASE("Containment", "[unary_ops][contains]", GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;

    for (size_t dim = dim_start; dim < dim_end; ++dim) {
        auto set = TestType::from_unit_box(dim);

        skip_if_unimplemented([&] {
            // Generate 10 random points within [-1, 1] for each dimension
            Eigen::MatrixXd samples =
                Eigen::MatrixXd::Random(10, dim); // Random returns values in [-1, 1]

            for (int i = 0; i < samples.rows(); ++i) {
                Eigen::VectorXd point = samples.row(i);

                bool contains = set.contains_point(point);
                REQUIRE(contains);
            }
        });
    }
}

TEMPLATE_LIST_TEST_CASE("Line intersection returns both unit box boundary points",
                        "[geosets][line]", GeoSetTypes) {
    const auto set = TestType::from_unit_box(2);

    skip_if_unimplemented([&] {
        const auto [lower, upper] =
            set.line_intersection(Eigen::Vector2d::Zero(), Eigen::Vector2d(1.0, 0.0));

        REQUIRE(lower.isApprox(Eigen::Vector2d(-1.0, 0.0)));
        REQUIRE(upper.isApprox(Eigen::Vector2d(1.0, 0.0)));
    });
}
