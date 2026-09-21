#include "catch2/catch_test_macros.hpp"
#include "geosets/sets/polygon.h"
#include <catch2/catch_all.hpp>
#include <geosets/parallel.h>
#include <geosets/sets/hpolytope.h>
#include <geosets/sets/interface.h>
#include <geosets/sets/vpolytope.h>

using GeoSetTypes = std::tuple<geosets::HPolytope, geosets::VPolytope, geosets::Polygon>;

TEMPLATE_LIST_TEST_CASE("Parallel execution with return type", "[geosets][parallel]", GeoSetTypes) {
    std::vector<TestType> polys;
    for (int i = 0; i < 10; ++i)
        polys.emplace_back(TestType::from_unit_box(2));

    geosets::set_omp_num_threads(8);
    auto centers = parallel_execute(polys, [](const TestType& p, size_t) { return p.center(); });
    REQUIRE(centers.size() == polys.size());
    for (const auto& c : centers) {
        REQUIRE(c.isApprox(Eigen::VectorXd::Zero(2).eval()));
    }
}

TEMPLATE_LIST_TEST_CASE("Parallel execution without return type", "[geosets][parallel]",
                        GeoSetTypes) {
    std::vector<TestType> polys;
    for (int i = 0; i < 10; ++i)
        polys.emplace_back(TestType::from_unit_box(2));

    geosets::set_omp_num_threads(8);
    parallel_execute_void(
        polys, [&](TestType& p, size_t) { p.translate_(Eigen::VectorXd::Ones(2).eval()); });
    for (size_t i = 0; i < polys.size(); ++i) {
        REQUIRE(polys[i].center().isApprox(Eigen::VectorXd::Ones(2).eval()));
    }
}

TEMPLATE_LIST_TEST_CASE("Parallel execution uses index argument", "[geosets][parallel]",
                        GeoSetTypes) {
    std::vector<TestType> polys;
    const size_t N = 20;
    for (size_t i = 0; i < N; ++i)
        polys.emplace_back(TestType::from_unit_box(2));

    geosets::set_omp_num_threads(4);
    // Translate each by (i, i)
    parallel_execute_void(polys, [](TestType& p, size_t i) {
        p.translate_(Eigen::VectorXd::Constant(2, static_cast<double>(i)).eval());
    });
    for (size_t i = 0; i < N; ++i) {
        REQUIRE(polys[i].center().isApprox(
            Eigen::VectorXd::Constant(2, static_cast<double>(i)).eval()));
    }
}

TEMPLATE_LIST_TEST_CASE("Parallel execution returns scalar", "[geosets][parallel]", GeoSetTypes) {
    std::vector<TestType> polys;
    const size_t N = 15;
    for (size_t i = 0; i < N; ++i)
        polys.emplace_back(TestType::from_unit_box(2));

    geosets::set_omp_num_threads(4);
    auto volumes = parallel_execute(polys, [](const TestType& p, size_t) { return p.volume(); });
    REQUIRE(volumes.size() == N);
    for (const auto& v : volumes) {
        REQUIRE(v == Catch::Approx(4.0));
    }
}

TEMPLATE_LIST_TEST_CASE("Parallel execution on empty vector", "[geosets][parallel]", GeoSetTypes) {
    std::vector<TestType> polys;
    geosets::set_omp_num_threads(4);
    auto centers = parallel_execute(polys, [](const TestType& p, size_t) { return p.center(); });
    REQUIRE(centers.empty());
    parallel_execute_void(polys, [](TestType& p, size_t) { (void)p; });
    REQUIRE(polys.empty());
}

TEMPLATE_LIST_TEST_CASE("Parallel execution single element", "[geosets][parallel]", GeoSetTypes) {
    std::vector<TestType> polys;
    polys.emplace_back(TestType::from_unit_box(2));

    geosets::set_omp_num_threads(4);
    auto dims = parallel_execute(polys, [](const TestType& p, size_t) { return p.dim(); });
    REQUIRE(dims.size() == 1);
    REQUIRE(dims[0] == 2);
}

TEMPLATE_LIST_TEST_CASE("Parallel execution single-threaded matches multi-threaded",
                        "[geosets][parallel]", GeoSetTypes) {
    std::vector<TestType> polys_a, polys_b;
    const size_t N = 25;
    for (size_t i = 0; i < N; ++i) {
        polys_a.emplace_back(TestType::from_unit_box(2));
        polys_b.emplace_back(TestType::from_unit_box(2));
    }

    geosets::set_omp_num_threads(1);
    auto centers_serial =
        parallel_execute(polys_a, [](const TestType& p, size_t) { return p.center(); });

    geosets::set_omp_num_threads(8);
    auto centers_parallel =
        parallel_execute(polys_b, [](const TestType& p, size_t) { return p.center(); });

    REQUIRE(centers_serial.size() == centers_parallel.size());
    for (size_t i = 0; i < centers_serial.size(); ++i) {
        REQUIRE(centers_serial[i].isApprox(centers_parallel[i]));
    }
}
