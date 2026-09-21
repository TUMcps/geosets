#include <geosets/sampling/gaussian_rejection.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

namespace geosets::sampling {

namespace {

void sample_diagonal_gaussian(const Eigen::VectorXd& loc, const Eigen::VectorXd& scale,
                              Eigen::VectorXd& sample) {
    for (Eigen::Index dim_idx = 0; dim_idx < loc.size(); ++dim_idx) {
        sample(dim_idx) = loc(dim_idx) + scale(dim_idx) * detail::sample_standard_normal();
    }
}

void throw_if_rejection_failed(const std::vector<std::uint8_t>& not_accepted,
                               std::size_t sampling_limit) {
    const auto failed = std::find(not_accepted.begin(), not_accepted.end(), 1);
    if (failed == not_accepted.end()) return;

    const std::size_t sample_idx = static_cast<std::size_t>(failed - not_accepted.begin());
    throw std::runtime_error("Gaussian rejection sampling failed for sample index " +
                             std::to_string(sample_idx) + " after " +
                             std::to_string(sampling_limit) + " resampling rounds.");
}

} // namespace

namespace detail {

RejectionResult gaussian_rejection_samples_with_mask(const GeometricSetInterface& set,
                                                     const Eigen::VectorXd& loc,
                                                     const Eigen::VectorXd& scale,
                                                     std::size_t samples_per_set,
                                                     std::size_t sampling_limit) {
    RejectionResult result{
        RowMajorMatrixXd(static_cast<Eigen::Index>(samples_per_set), loc.size()),
        std::vector<std::uint8_t>(samples_per_set, 1),
    };
    std::size_t remaining = samples_per_set;
    Eigen::VectorXd candidate(loc.size());

    const auto try_accept = [&](std::size_t sample_idx) {
        sample_diagonal_gaussian(loc, scale, candidate);
        if (!set.contains_point(candidate, 1e-9, false)) return;

        result.samples.row(static_cast<Eigen::Index>(sample_idx)) = candidate.transpose();
        result.not_accepted[sample_idx] = 0;
        --remaining;
    };

    // One initial draw is always performed before sampling_limit resampling rounds.
    for (std::size_t sample_idx = 0; sample_idx < samples_per_set; ++sample_idx) {
        try_accept(sample_idx);
    }

    for (std::size_t round = 0; round < sampling_limit && remaining != 0; ++round) {
        for (std::size_t sample_idx = 0; sample_idx < samples_per_set; ++sample_idx) {
            if (result.not_accepted[sample_idx] == 0) continue;
            try_accept(sample_idx);
        }
    }

    return result;
}

} // namespace detail

RowMajorMatrixXd gaussian_rejection_samples(const GeometricSetInterface& set,
                                            const Eigen::VectorXd& loc,
                                            const Eigen::VectorXd& scale, std::size_t n_samples,
                                            std::size_t rejection_limit, bool validate) {
    detail::validate_sample_inputs(set, loc, scale, validate);
    detail::RejectionResult result =
        detail::gaussian_rejection_samples_with_mask(set, loc, scale, n_samples, rejection_limit);
    throw_if_rejection_failed(result.not_accepted, rejection_limit);
    return std::move(result.samples);
}

} // namespace geosets::sampling
