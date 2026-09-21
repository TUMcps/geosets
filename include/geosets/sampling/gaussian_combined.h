#pragma once

#include <geosets/sampling/common.h>

#include <cstddef>

namespace geosets::sampling {

RowMajorMatrixXd gaussian_rejection_rdhr_fallback_samples(
    const GeometricSetInterface& set, const Eigen::VectorXd& loc, const Eigen::VectorXd& scale,
    std::size_t n_samples = 1, std::size_t rejection_limit = 100, bool validate = true);

} // namespace geosets::sampling
