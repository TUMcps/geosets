#include "../definitions.h"

#include <geosets/core.h>
#include <geosets/sampling/common.h>
#include <geosets/sampling/gaussian_combined.h>
#include <geosets/sampling/gaussian_rdhr.h>
#include <geosets/sampling/gaussian_rejection.h>

#include <catch2/catch_all.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using geosets::sampling::RowMajorMatrixXd;

namespace {

// Uniform adaptor over the three samplers so a single common test can exercise all of them.
// Generic lambdas bind each function's default rejection_limit / validate.
using Sampler =
    std::function<RowMajorMatrixXd(const geosets::GeometricSetInterface&, const Eigen::VectorXd&,
                                   const Eigen::VectorXd&, std::size_t)>;

const std::vector<std::pair<std::string, Sampler>> kSamplers = {
    {"rejection",
     [](const geosets::GeometricSetInterface& s, const Eigen::VectorXd& loc,
        const Eigen::VectorXd& scale,
        std::size_t n) { return geosets::sampling::gaussian_rejection_samples(s, loc, scale, n); }},
    {"rdhr",
     [](const geosets::GeometricSetInterface& s, const Eigen::VectorXd& loc,
        const Eigen::VectorXd& scale,
        std::size_t n) { return geosets::sampling::gaussian_rdhr_samples(s, loc, scale, n); }},
    {"combined",
     [](const geosets::GeometricSetInterface& s, const Eigen::VectorXd& loc,
        const Eigen::VectorXd& scale, std::size_t n) {
         return geosets::sampling::gaussian_rejection_rdhr_fallback_samples(s, loc, scale, n);
     }},
};

void require_contained(const RowMajorMatrixXd& samples, const geosets::GeometricSetInterface& set,
                       std::size_t n, double tol = 1e-6) {
    REQUIRE(samples.rows() == static_cast<Eigen::Index>(n));
    REQUIRE(samples.cols() == static_cast<Eigen::Index>(set.dim()));
    for (Eigen::Index row = 0; row < samples.rows(); ++row) {
        // MCMC samples can land essentially on the boundary; a loosened tol avoids false fails.
        REQUIRE(set.contains_point(samples.row(row).transpose(), tol));
    }
}

} // namespace

// Every returned sample must actually lie inside the set (real correctness for the MCMC paths,
// not just a shape check), and n == 0 must yield an empty (0 x dim) matrix.
TEMPLATE_LIST_TEST_CASE("Gaussian sampling returns contained samples of correct shape",
                        "[sampling][gaussian][containment]", GeoSetTypes) {
    auto [dim_start, dim_end] = DefaultDimensions<TestType>::range;
    for (std::size_t dim = dim_start; dim < dim_end; ++dim) {
        auto set = TestType::from_unit_box(dim);
        const Eigen::VectorXd loc = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(dim));
        const Eigen::VectorXd scale =
            Eigen::VectorXd::Constant(static_cast<Eigen::Index>(dim), 0.5);

        for (const auto& entry : kSamplers) {
            const std::string& name = entry.first;
            const Sampler& sample = entry.second;
            INFO(name << " dim=" << dim);
            skip_if_unimplemented([&] {
                geosets::random_seed(17);
                RowMajorMatrixXd samples = sample(set, loc, scale, 16);
                require_contained(samples, set, 16);

                RowMajorMatrixXd none =
                    sample(set, loc, Eigen::VectorXd::Ones(static_cast<Eigen::Index>(dim)), 0);
                REQUIRE(none.rows() == 0);
                REQUIRE(none.cols() == static_cast<Eigen::Index>(dim));
            });
        }
    }
}

// The one distribution-level test. The unit box and a loc=0 Gaussian are both symmetric under
// x -> -x, and RDHR starts from the (symmetric) center, so every sample's marginal is mean-zero
// for ANY walk length -- only sampling noise remains. This sidesteps RDHR's weak mixing
// (walk_length = dim^3 = 8 at dim 2) which would make an off-center mean check flaky.
TEMPLATE_LIST_TEST_CASE("Gaussian sampling of a symmetric set is empirically centered",
                        "[sampling][gaussian][statistical]", GeoSetTypes) {
    const std::size_t dim = 2;
    auto set = TestType::from_unit_box(dim);
    const Eigen::VectorXd loc = Eigen::VectorXd::Zero(dim);
    const Eigen::VectorXd scale = Eigen::VectorXd::Constant(dim, 0.6);
    const std::size_t n = 256;

    for (const auto& entry : kSamplers) {
        const std::string& name = entry.first;
        const Sampler& sample = entry.second;
        skip_if_unimplemented([&] {
            geosets::random_seed(101);
            RowMajorMatrixXd samples = sample(set, loc, scale, n);
            Eigen::RowVectorXd mean = samples.colwise().mean();
            for (Eigen::Index j = 0; j < mean.size(); ++j) {
                INFO(name << " coord " << j << " mean=" << mean(j));
                REQUIRE(std::abs(mean(j)) < 0.12); // ~4 * standard error
            }
        });
    }
}

// The only source of nondeterminism is the global RNG: same seed => identical output.
TEMPLATE_LIST_TEST_CASE("Gaussian sampling is reproducible under a fixed seed",
                        "[sampling][gaussian][determinism]", GeoSetTypes) {
    const std::size_t dim = 2;
    auto set = TestType::from_unit_box(dim);
    const Eigen::VectorXd loc = Eigen::VectorXd::Zero(dim);
    const Eigen::VectorXd scale = Eigen::VectorXd::Constant(dim, 0.5);

    for (const auto& entry : kSamplers) {
        const std::string& name = entry.first;
        const Sampler& sample = entry.second;
        INFO(name);
        skip_if_unimplemented([&] {
            geosets::random_seed(7);
            RowMajorMatrixXd a = sample(set, loc, scale, 8);
            geosets::random_seed(7);
            RowMajorMatrixXd b = sample(set, loc, scale, 8);
            REQUIRE(a.isApprox(b));
        });
    }
}

// When acceptance is easy no fallback fires, so the combined sampler must consume the RNG stream
// exactly like pure rejection and produce a byte-for-byte identical result.
TEMPLATE_LIST_TEST_CASE("Combined sampling equals rejection when no fallback is needed",
                        "[sampling][gaussian][combined]", GeoSetTypes) {
    const std::size_t dim = 2;
    auto set = TestType::from_unit_box(dim);
    const Eigen::VectorXd loc = Eigen::VectorXd::Zero(dim);
    const Eigen::VectorXd scale = Eigen::VectorXd::Constant(dim, 0.5);
    const std::size_t n = 8;
    const std::size_t limit = 10000;

    skip_if_unimplemented([&] {
        geosets::random_seed(31);
        RowMajorMatrixXd rejection =
            geosets::sampling::gaussian_rejection_samples(set, loc, scale, n, limit);
        geosets::random_seed(31);
        RowMajorMatrixXd combined =
            geosets::sampling::gaussian_rejection_rdhr_fallback_samples(set, loc, scale, n, limit);
        REQUIRE(combined.isApprox(rejection, 1e-12));
    });
}

// Differential behaviour that justifies three functions: with the Gaussian mass parked far outside
// the set and a tiny budget, exact rejection gives up (throws) while RDHR and the combined
// fallback still deliver contained samples via the interior-point (center) fallback.
TEMPLATE_LIST_TEST_CASE("Rejection exhausts its budget where RDHR and combined survive",
                        "[sampling][gaussian][fallback]", GeoSetTypes) {
    const std::size_t dim = 2;
    auto set = TestType::from_unit_box(dim);
    const Eigen::VectorXd loc = Eigen::VectorXd::Constant(dim, 5.0); // outside [-1, 1]^2
    const Eigen::VectorXd scale = Eigen::VectorXd::Constant(dim, 0.1);
    const std::size_t n = 8;

    skip_if_unimplemented([&] {
        geosets::random_seed(41);
        REQUIRE_THROWS_AS(geosets::sampling::gaussian_rejection_samples(set, loc, scale, n, 1),
                          std::runtime_error);

        RowMajorMatrixXd rdhr = geosets::sampling::gaussian_rdhr_samples(set, loc, scale, n);
        require_contained(rdhr, set, n);

        RowMajorMatrixXd combined =
            geosets::sampling::gaussian_rejection_rdhr_fallback_samples(set, loc, scale, n, 1);
        require_contained(combined, set, n);
    });
}

// Shared input validation (all three route through validate_sample_inputs).
TEMPLATE_LIST_TEST_CASE("Gaussian sampling rejects invalid inputs",
                        "[sampling][gaussian][validation]", GeoSetTypes) {
    const std::size_t dim = 2;
    const Eigen::VectorXd loc = Eigen::VectorXd::Zero(dim);
    const Eigen::VectorXd scale = Eigen::VectorXd::Ones(dim);

    for (const auto& entry : kSamplers) {
        const std::string& name = entry.first;
        const Sampler& sample = entry.second;
        INFO(name);

        auto box = TestType::from_unit_box(dim);

        // Infinite / empty sets are unsupported.
        REQUIRE_THROWS(sample(TestType::make_infinite(dim), loc, scale, 4));
        REQUIRE_THROWS(sample(TestType::make_empty(dim), loc, scale, 4));

        // Non-finite loc and non-positive scale are rejected.
        Eigen::VectorXd bad_loc = loc;
        bad_loc(0) = std::numeric_limits<double>::infinity();
        REQUIRE_THROWS(sample(box, bad_loc, scale, 4));

        Eigen::VectorXd bad_scale = scale;
        bad_scale(0) = 0.0;
        REQUIRE_THROWS(sample(box, loc, bad_scale, 4));

        // loc / scale dimension mismatch is rejected even structurally.
        REQUIRE_THROWS(sample(box, Eigen::VectorXd::Zero(dim + 1), scale, 4));
        REQUIRE_THROWS(sample(box, loc, Eigen::VectorXd::Ones(dim + 1), 4));
    }
}
