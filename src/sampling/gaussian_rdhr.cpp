#include <geosets/sampling/gaussian_rdhr.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <boost/math/distributions/normal.hpp>

namespace geosets::sampling {

namespace {

Eigen::VectorXd random_unit_direction(Eigen::Index dim) {
    Eigen::VectorXd direction(dim);
    double norm = 0.0;
    do {
        for (Eigen::Index i = 0; i < dim; ++i) {
            direction(i) = detail::sample_standard_normal();
        }
        norm = direction.norm();
    } while (norm <= 1e-12);
    return direction / norm;
}

// Note: Two options here for random walk step:
//  1. propose with line-segment-truncated gaussian and no metropolis filter.
//  2. propose with uniform distribution over line segment and reject/accept with metropolis filter.
// Currently we use 2., since it better handles numerical instabilities with very low probability
// mass in the set.
[[maybe_unused]] double sample_truncated_gaussian_1d(double lower, double upper, double mean,
                                                     double sigma, double eps = 1e-5) {
    if (!std::isfinite(mean)) {
        throw std::runtime_error("Truncated Gaussian mean must be finite.");
    }
    if (!(sigma > 0.0) || !std::isfinite(sigma)) {
        throw std::runtime_error("Truncated Gaussian sigma must be strictly positive and finite.");
    }

    const double shrunk_lower = lower + eps;
    const double shrunk_upper = upper - eps;
    if (!(shrunk_lower < shrunk_upper)) {
        throw std::runtime_error("Truncated Gaussian interval is too narrow after epsilon shrink.");
    }

    const boost::math::normal_distribution<double> distribution(mean, sigma);
    const double cdf_lower = boost::math::cdf(distribution, lower);
    const double cdf_upper = boost::math::cdf(distribution, upper);
    if (!std::isfinite(cdf_lower) || !std::isfinite(cdf_upper) || !(cdf_lower < cdf_upper)) {
        throw std::runtime_error("Truncated Gaussian CDF range collapsed.");
    }

    const double probability =
        cdf_lower + (cdf_upper - cdf_lower) * detail::sample_unit_uniform_open();
    const double sample = boost::math::quantile(distribution, probability);
    if (!std::isfinite(sample)) {
        throw std::runtime_error("Truncated Gaussian quantile produced a non-finite sample.");
    }
    return std::clamp(sample, shrunk_lower, shrunk_upper);
}

[[maybe_unused]] double sample_gaussian_line_parameter(const std::pair<double, double>& bounds,
                                                       const Eigen::VectorXd& direction,
                                                       const Eigen::VectorXd& point,
                                                       const Eigen::VectorXd& loc,
                                                       const Eigen::VectorXd& scale) {
    const double mean = direction.dot(loc - point);
    const double sigma = (scale.array() * direction.array()).matrix().norm();
    return sample_truncated_gaussian_1d(bounds.first, bounds.second, mean, sigma);
}

[[maybe_unused]] double sample_uniform_line_parameter(const std::pair<double, double>& bounds) {
    return bounds.first + (bounds.second - bounds.first) * detail::sample_unit_uniform_open();
}

double gaussian_log_kernel(const Eigen::VectorXd& point, const Eigen::VectorXd& loc,
                           const Eigen::VectorXd& scale) {
    return -0.5 * ((point - loc).array() / scale.array()).matrix().squaredNorm();
}

[[maybe_unused]] bool metropolis_accept(const Eigen::VectorXd& current,
                                        const Eigen::VectorXd& proposal, const Eigen::VectorXd& loc,
                                        const Eigen::VectorXd& scale) {
    const double log_accept_ratio = std::min(
        gaussian_log_kernel(proposal, loc, scale) - gaussian_log_kernel(current, loc, scale), 0.0);
    return log_accept_ratio >= std::log(detail::sample_unit_uniform_open());
}

void rdhr_step(const GeometricSetInterface& set, Eigen::VectorXd& current,
               const Eigen::VectorXd& loc, const Eigen::VectorXd& scale) {
    const Eigen::VectorXd direction = random_unit_direction(current.size());

    const auto bounds = set.line_intersection_bounds(current, direction, false);
    if (!std::isfinite(bounds.first) || !std::isfinite(bounds.second) ||
        !(bounds.first < bounds.second)) {
        throw std::runtime_error(
            "RDHR step failed to find valid line intersection bounds for the current point.");
    }

    // const double line_parameter =
    //     sample_gaussian_line_parameter(bounds, direction, current, loc, scale);
    const double line_parameter = sample_uniform_line_parameter(bounds);
    Eigen::VectorXd proposal = current + line_parameter * direction;
    if (metropolis_accept(current, proposal, loc, scale)) {
        current.swap(proposal);
    }
}

} // namespace

namespace detail {

Eigen::VectorXd gaussian_rdhr_initial_point(const GeometricSetInterface& set,
                                            const Eigen::VectorXd& loc) {
    Eigen::VectorXd initial = loc;
    if (!set.contains_point(initial, 1e-9, false)) initial = set.center(false);
    return initial;
}

Eigen::VectorXd gaussian_rdhr_from_initial(const GeometricSetInterface& set,
                                           const Eigen::VectorXd& loc, const Eigen::VectorXd& scale,
                                           const Eigen::VectorXd& initial) {
    Eigen::VectorXd current = initial;
    size_t dim = set.dim();
    const std::size_t walk_length = dim * dim * dim;

    for (std::size_t step = 0; step < walk_length; ++step) {
        rdhr_step(set, current, loc, scale);
    }

    return current;
}

} // namespace detail

RowMajorMatrixXd gaussian_rdhr_samples(const GeometricSetInterface& set, const Eigen::VectorXd& loc,
                                       const Eigen::VectorXd& scale, std::size_t n_samples,
                                       bool validate) {
    detail::validate_sample_inputs(set, loc, scale, validate);

    RowMajorMatrixXd samples(static_cast<Eigen::Index>(n_samples), loc.size());
    if (n_samples == 0) return samples;
    const Eigen::VectorXd initial = detail::gaussian_rdhr_initial_point(set, loc);
    for (std::size_t sample_idx = 0; sample_idx < n_samples; ++sample_idx) {
        samples.row(static_cast<Eigen::Index>(sample_idx)) =
            detail::gaussian_rdhr_from_initial(set, loc, scale, initial).transpose();
    }
    return samples;
}

} // namespace geosets::sampling
