#include <geosets/sampling/gaussian_combined.h>

#include <geosets/sampling/gaussian_rdhr.h>
#include <geosets/sampling/gaussian_rejection.h>

#include <utility>

namespace geosets::sampling {

RowMajorMatrixXd gaussian_rejection_rdhr_fallback_samples(
    const GeometricSetInterface& set, const Eigen::VectorXd& loc, const Eigen::VectorXd& scale,
    std::size_t n_samples, std::size_t rejection_limit, bool validate) {
    detail::validate_sample_inputs(set, loc, scale, validate);
    detail::RejectionResult result =
        detail::gaussian_rejection_samples_with_mask(set, loc, scale, n_samples, rejection_limit);

    Eigen::VectorXd initial;
    for (std::size_t sample_idx = 0; sample_idx < n_samples; ++sample_idx) {
        if (result.not_accepted[sample_idx] == 0) continue;
        if (initial.size() == 0) {
            initial = detail::gaussian_rdhr_initial_point(set, loc);
        }
        result.samples.row(static_cast<Eigen::Index>(sample_idx)) =
            detail::gaussian_rdhr_from_initial(set, loc, scale, initial).transpose();
    }
    return std::move(result.samples);
}

} // namespace geosets::sampling
