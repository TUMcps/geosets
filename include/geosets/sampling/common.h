#pragma once

#include <geosets/sets/interface.h>

#include <Eigen/Dense>

#include <cstddef>

namespace geosets::sampling {

using RowMajorMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

namespace detail {

void validate_sample_inputs(const GeometricSetInterface& set, const Eigen::VectorXd& loc,
                            const Eigen::VectorXd& scale, bool validate);

double sample_unit_uniform_open();
double sample_standard_normal();

} // namespace detail

} // namespace geosets::sampling
