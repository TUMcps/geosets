#pragma once

#include <geosets/sampling/common.h>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace geosets::sampling {

namespace detail {

struct RejectionResult {
    RowMajorMatrixXd samples;
    std::vector<std::uint8_t> not_accepted;
};

RejectionResult gaussian_rejection_samples_with_mask(const GeometricSetInterface& set,
                                                     const Eigen::VectorXd& loc,
                                                     const Eigen::VectorXd& scale,
                                                     std::size_t samples_per_set,
                                                     std::size_t sampling_limit);

} // namespace detail

inline constexpr size_t DEFAULT_REJECTION_LIMIT = 10000;

RowMajorMatrixXd gaussian_rejection_samples(const GeometricSetInterface& set,
                                            const Eigen::VectorXd& loc,
                                            const Eigen::VectorXd& scale, std::size_t n_samples = 1,
                                            std::size_t rejection_limit = DEFAULT_REJECTION_LIMIT,
                                            bool validate = true);

} // namespace geosets::sampling
