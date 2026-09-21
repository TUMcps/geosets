#pragma once

#include <Eigen/Dense>
#include <Highs.h>

namespace geosets::highs_wrapper {

// If check_scaling is true, reject inputs with a constraint row whose coefficients all fall at or
// below HiGHS's small_matrix_value tolerance (such a row would be silently culled, changing the
// problem). Off by default: the caller opts in when it wants that safety over raw throughput.
HighsModelStatus solve_lp(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                          const Eigen::VectorXd& objective, Eigen::VectorXd* solution = nullptr,
                          bool check_scaling = false);

inline bool is_unbounded_status(HighsModelStatus status) {
    return status == HighsModelStatus::kUnbounded ||
           status == HighsModelStatus::kUnboundedOrInfeasible;
}

inline bool is_infeasible_status(HighsModelStatus status) {
    return status == HighsModelStatus::kInfeasible ||
           status == HighsModelStatus::kUnboundedOrInfeasible;
}

} // namespace geosets::highs_wrapper
