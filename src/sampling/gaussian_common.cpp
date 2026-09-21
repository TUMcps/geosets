#include <geosets/sampling/common.h>

#include <geosets/exceptions.h>
#include <geosets/parallel.h>

#include <limits>
#include <random>
#include <stdexcept>

namespace geosets::sampling::detail {

namespace {

void validate_distribution_parameters(const Eigen::VectorXd& loc, const Eigen::VectorXd& scale) {
    if (!loc.array().isFinite().all() || !scale.array().isFinite().all()) {
        throw std::invalid_argument("Gaussian loc and scale must contain only finite values.");
    }
    if (!(scale.array() > 0.0).all()) {
        throw std::invalid_argument("Gaussian scale entries must be strictly positive.");
    }
}

void validate_set(const GeometricSetInterface& set) {
    if (set.is_infinite()) {
        throw InfiniteSetException("Gaussian constrained sampling requires a bounded set.");
    }
    if (set.is_empty()) {
        throw EmptySetException("Gaussian constrained sampling cannot use an empty set.");
    }
    if (!set.bounded()) {
        throw std::invalid_argument("Gaussian constrained sampling requires a bounded set.");
    }
    if (set.degenerate()) {
        throw DegenerateSetException();
    }
}

} // namespace

void validate_sample_inputs(const GeometricSetInterface& set, const Eigen::VectorXd& loc,
                            const Eigen::VectorXd& scale, bool validate) {
    // These structural checks protect the unchecked sampling helpers.
    const Eigen::Index dim = static_cast<Eigen::Index>(set.dim());
    if (dim == 0) {
        throw std::invalid_argument("Gaussian sampling requires a positive set dimension.");
    }
    if (loc.size() != dim || scale.size() != dim) {
        throw std::invalid_argument("Gaussian loc and scale must match the set dimension.");
    }

    // Trusted callers can skip the more expensive semantic checks.
    if (!validate) return;
    validate_distribution_parameters(loc, scale);
    validate_set(set);
}

double sample_unit_uniform_open() {
    // Open interval: lower bound off zero keeps callers' std::log finite (never log(0)).
    std::uniform_real_distribution<double> dist(std::numeric_limits<double>::min(), 1.0);
    return dist(thread_rng());
}

double sample_standard_normal() {
    // Fresh distribution each call: normal_distribution caches a spare Box-Muller value, which
    // would leak across a reseed and break reproducibility.
    std::normal_distribution<double> dist;
    return dist(thread_rng());
}

} // namespace geosets::sampling::detail
