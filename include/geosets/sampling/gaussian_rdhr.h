#pragma once

#include <geosets/sampling/common.h>

#include <cstddef>

namespace geosets::sampling {

namespace detail {

Eigen::VectorXd gaussian_rdhr_initial_point(const GeometricSetInterface& set,
                                            const Eigen::VectorXd& loc);
Eigen::VectorXd gaussian_rdhr_from_initial(const GeometricSetInterface& set,
                                           const Eigen::VectorXd& loc, const Eigen::VectorXd& scale,
                                           const Eigen::VectorXd& initial);

} // namespace detail

RowMajorMatrixXd gaussian_rdhr_samples(const GeometricSetInterface& set, const Eigen::VectorXd& loc,
                                       const Eigen::VectorXd& scale, std::size_t n_samples = 1,
                                       bool validate = true);

} // namespace geosets::sampling
