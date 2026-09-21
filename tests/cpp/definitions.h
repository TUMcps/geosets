#include "catch2/catch_test_macros.hpp"
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <catch2/catch_all.hpp>
#include <geosets/sets/hpolytope.h>
#include <geosets/sets/interval.h>
#include <geosets/sets/polygon.h>
#include <geosets/sets/vpolytope.h>
#include <geosets/sets/zonotope.h>

using GeoSetTypes = std::tuple<geosets::HPolytope, geosets::VPolytope, geosets::Polygon,
                               geosets::Interval, geosets::Zonotope>;

constexpr double tol_ = 1e-9;

// -- Dimension ranges for tests --
template <typename T> struct DefaultDimensions {
    static constexpr std::pair<std::size_t, std::size_t> range{2, 5}; // [begin, end)
};

template <> struct DefaultDimensions<geosets::Polygon> {
    static constexpr std::pair<std::size_t, std::size_t> range{2, 3};
};

template <> struct DefaultDimensions<geosets::Zonotope> {
    static constexpr std::pair<std::size_t, std::size_t> range{2, 5};
};

// Macro for detecting not implemented features
template <typename F> void skip_if_unimplemented(F&& fn) {
    try {
        fn();
    } catch (const geosets::NotImplemented&) {
        SKIP("Not implemented");
    }
}

#define REQUIRE_OR_SKIP(expr)                                                                      \
    do {                                                                                           \
        try {                                                                                      \
            auto _result = (expr);                                                                 \
            REQUIRE(_result);                                                                      \
        } catch (const geosets::NotImplemented&) {                                                 \
            SKIP("Not implemented");                                                               \
        }                                                                                          \
    } while (false)

#define REQUIRE_FALSE_OR_SKIP(expr)                                                                \
    do {                                                                                           \
        try {                                                                                      \
            auto _result = (expr);                                                                 \
            REQUIRE_FALSE(_result);                                                                \
        } catch (const geosets::NotImplemented&) {                                                 \
            SKIP("Not implemented");                                                               \
        }                                                                                          \
    } while (false)

#define REQUIRE_THROWS_AS_OR_SKIP(expr, exception_type)                                            \
    do {                                                                                           \
        try {                                                                                      \
            (expr);                                                                                \
            FAIL("Expected exception of type " #exception_type);                                   \
        } catch (const exception_type&) {                                                          \
        } catch (const geosets::NotImplemented&) {                                                 \
            SKIP("Not implemented");                                                               \
        }                                                                                          \
    } while (false)
