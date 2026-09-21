#include "geosets/structs.h"
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <Highs.h>
#include <cstddef>
#include <cstdio>
#include <geosets/conversion/to_vpolytope.h>
#include <geosets/exceptions.h>
#include <geosets/highs_wrapper.h>
#include <geosets/infinite.h>
#include <geosets/linalg_utils.h>
#include <geosets/sets/hpolytope.h>

namespace {

/**
 * @brief Checks if the point is on the active side of at least one inequality constraint (i.e., if
 * it is not strictly in the interior). This is used to determine if a point is on the boundary of
 * the HPolytope.
 */
bool has_active_inequality(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                           const Eigen::VectorXd& point, double tol) {
    Eigen::VectorXd residual = A * point - b;
    return (residual.array().abs() <= tol).any();
}

} // namespace

namespace geosets {
// --- Static constructors ---

HPolytope HPolytope::from_unit_box(size_t dim) {
    auto A =
        linalg::vstack(Eigen::MatrixXd::Identity(dim, dim), -Eigen::MatrixXd::Identity(dim, dim));
    auto b = Eigen::VectorXd::Ones(2 * dim);
    return HPolytope(std::move(A), std::move(b));
}

HPolytope HPolytope::from_random(size_t dim, size_t n_constraints) {
    auto box_oly = HPolytope::from_unit_box(dim);

    auto random_A = Eigen::MatrixXd::Random(n_constraints, dim).eval();

    // Normalize A
    for (size_t i = 0; i < n_constraints; ++i) {
        random_A.row(i).normalize();
    }

    // Uniform in (-0.8, 0.8)
    auto interior_point = Eigen::VectorXd::Random(dim) * 0.8;
    // Random positive offsets in (0.1, 0.3) to ensure interior_point satisfies all constraints
    Eigen::VectorXd offsets = (Eigen::VectorXd::Random(n_constraints).array() + 1.0) * 0.1 + 0.1;
    auto random_b = random_A * interior_point + offsets;

    // Concatenate to A
    auto A = linalg::vstack(box_oly.A, random_A);
    auto b = linalg::vstack(box_oly.b, random_b);

    return HPolytope(std::move(A), std::move(b));
}

// --- Properties ---

std::string HPolytope::str() const {
    std::stringstream ss;
    Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n", "[", "]");
    ss << "A:\n" << this->A.format(fmt) << "\n";
    ss << "b:\n" << this->b.format(fmt);
    return "HPolytope\n" + ss.str();
}

bool HPolytope::degenerate_impl() const {
    Eigen::VectorXd center;
    try {
        center = this->center();
    } catch (const FailedOptimizationProblem&) {
        return true;
    }
    return has_active_inequality(this->A, this->b, center, this->atol_);
}

bool HPolytope::is_empty_impl() const {
    const Eigen::VectorXd objective = Eigen::VectorXd::Zero(this->dim());
    const HighsModelStatus status = highs_wrapper::solve_lp(this->A, this->b, objective);
    return highs_wrapper::is_infeasible_status(status);
}

bool HPolytope::bounded_impl() const {
    const size_t n = this->dim();
    for (size_t i = 0; i < n; ++i) {
        Eigen::VectorXd axis = Eigen::VectorXd::Zero(n);
        axis(i) = 1.0;

        const HighsModelStatus max_status = highs_wrapper::solve_lp(this->A, this->b, axis);
        if (highs_wrapper::is_unbounded_status(max_status)) {
            return false;
        }

        const HighsModelStatus min_status = highs_wrapper::solve_lp(this->A, this->b, -axis);
        if (highs_wrapper::is_unbounded_status(min_status)) {
            return false;
        }
    }

    return true;
}

// --- Rest ---

bool HPolytope::contains_point_impl(const Eigen::VectorXd& point, double tol) const {
    const Eigen::VectorXd residual = this->A * point - this->b;
    return residual.maxCoeff() <= tol;
}

Eigen::MatrixXd HPolytope::to_vertices_impl() const {
    return geosets::conversion::hpolytope_to_vpolytope(*this).vertices;
}

Eigen::VectorXd HPolytope::center_impl() const {
    if (this->cached_center_) return *this->cached_center_;

    // Chebyshev center: maximize r such that A*x + r*||a_i|| <= b_i
    // Formulated as LP with variables [x; r]
    size_t n = this->dim();
    size_t m = this->n_constraints();

    // Augment A with column of row norms for the radius variable
    Eigen::MatrixXd A_extended(m, n + 1);
    A_extended.leftCols(n) = this->A;

    for (size_t i = 0; i < m; ++i) {
        A_extended(i, n) = this->A.row(i).norm();
    }

    Eigen::VectorXd objective = Eigen::VectorXd::Zero(n + 1);
    objective(n) = 1.0;

    Eigen::VectorXd solution;
    const HighsModelStatus status =
        highs_wrapper::solve_lp(A_extended, this->b, objective, &solution, true);

    if (status != HighsModelStatus::kOptimal) {
        throw FailedOptimizationProblem();
    }

    this->cached_center_ = solution.head(n);
    return *this->cached_center_;
}

geosets::Support HPolytope::support_function_impl(const Eigen::VectorXd& direction) const {
    // For an H-polytope, the support function in direction d is:
    // max{ d^T x | Ax <= b }
    // Note: highs_wrapper::solve_lp internally negates the objective, so it maximizes
    Eigen::VectorXd solution;
    const HighsModelStatus status = highs_wrapper::solve_lp(this->A, this->b, direction, &solution);

    if (highs_wrapper::is_unbounded_status(status)) {
        // If unbounded, the support function value is infinite
        return {std::numeric_limits<double>::infinity(), Eigen::VectorXd()};
    }

    if (status != HighsModelStatus::kOptimal) {
        throw FailedOptimizationProblem();
    }

    double support_value = direction.dot(solution);
    return {support_value, solution};
}

geosets::Hyperplane HPolytope::tangent_hyperplane_impl(const Eigen::VectorXd& direction) const {
    Eigen::VectorXd c = this->center_impl();

    double norm_d = direction.norm();
    if (norm_d < this->atol_) {
        throw std::runtime_error(
            "Direction vector is too close to zero for tangent hyperplane computation.");
    }

    Eigen::VectorXd Ad = this->A * direction;
    Eigen::VectorXd slack = this->b - (this->A * c);

    double t_min = std::numeric_limits<double>::infinity();
    size_t min_index = 0;
    for (int i = 0; i < Ad.size(); ++i) {
        if (Ad(i) > this->atol_) {
            // t_min = std::min(t_min, slack(i) / Ad(i));
            if (slack(i) / Ad(i) < t_min) {
                t_min = slack(i) / Ad(i);
                min_index = static_cast<size_t>(i);
            }
        }
    }
    return {this->A.row(min_index), this->b(min_index)};
}

double HPolytope::volume_impl() const {
    return conversion::hpolytope_to_vpolytope(*this).volume();
}

Eigen::VectorXd HPolytope::boundary_point_impl(const Eigen::VectorXd& direction) const {
    Eigen::VectorXd c = this->center_impl();
    Eigen::VectorXd Ad = this->A * direction;
    Eigen::VectorXd slack = this->b - this->A * c;
    double t_min = std::numeric_limits<double>::infinity();

    for (int i = 0; i < Ad.size(); ++i) {
        if (Ad(i) > 0) {
            t_min = std::min(t_min, slack(i) / Ad(i));
        }
    }

    if (std::isinf(t_min)) return c;
    return c + t_min * direction;
}

std::pair<double, double>
HPolytope::line_intersection_bounds_impl(const Eigen::VectorXd& point,
                                         const Eigen::VectorXd& direction) const {
    const Eigen::VectorXd line_coefficients = this->A * direction;
    const Eigen::VectorXd slack = this->b - this->A * point;

    double lower = -std::numeric_limits<double>::infinity();
    double upper = std::numeric_limits<double>::infinity();

    for (Eigen::Index i = 0; i < line_coefficients.size(); ++i) {
        const double coefficient = line_coefficients(i);
        if (std::abs(coefficient) <= this->atol_) {
            if (slack(i) < -this->atol_) {
                throw std::runtime_error("The line does not intersect the HPolytope.");
            }
            continue;
        }

        const double boundary_parameter = slack(i) / coefficient;
        if (coefficient > 0.0) {
            upper = std::min(upper, boundary_parameter);
        } else {
            lower = std::max(lower, boundary_parameter);
        }
    }

    if (!std::isfinite(lower) || !std::isfinite(upper)) {
        throw std::runtime_error(
            "The line intersection is unbounded and has fewer than two boundary points.");
    }
    if (lower > upper + this->atol_) {
        throw std::runtime_error("The line does not intersect the HPolytope.");
    }

    return {lower, upper};
}

void HPolytope::translate_impl_(const Eigen::VectorXd& vector) {
    this->b = this->b + this->A * vector;
    if (this->cached_center_) *this->cached_center_ += vector;
}

void HPolytope::linear_transform_impl_(const Eigen::MatrixXd& matrix) {
    this->cached_center_.reset();
    const size_t m = matrix.rows();
    const size_t n = matrix.cols();
    const auto rank = linalg::rank(matrix, this->atol_);

    // Case 1: Full rank square matrix - use inverse
    if (m == n && rank == n) {
        this->A = this->A * matrix.inverse();
        return;
    }

    if (m > n) {
        throw std::runtime_error(
            "HPolytope linear transformation with dimension increase is not supported.");
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Eigen::MatrixXd U = svd.matrixU();
    const Eigen::MatrixXd V = svd.matrixV();
    const Eigen::VectorXd singular_values = svd.singularValues();

    Eigen::MatrixXd AV = this->A * V;

    const size_t r = static_cast<size_t>(rank);
    const size_t n_aux = n - r;

    // Variables are [y'; z], where y' = U^T y and z are eliminated auxiliary dimensions.
    // For i < r: x'_i = y'_i / sigma_i. For i >= r: x'_i = z_{i-r}.
    Eigen::MatrixXd A_new = Eigen::MatrixXd::Zero(this->A.rows(), m + n_aux);

    for (size_t i = 0; i < r; ++i) {
        A_new.col(i) = AV.col(static_cast<Eigen::Index>(i)) / singular_values(i);
    }
    if (n_aux > 0) {
        A_new.block(0, m, this->A.rows(), n_aux) = AV.rightCols(n_aux);
    }

    Eigen::VectorXd b_new = this->b;

    // For zero singular values, enforce y'_i = 0 via paired inequalities.
    if (r < m) {
        Eigen::MatrixXd A_eq = Eigen::MatrixXd::Zero(2 * (m - r), m + n_aux);
        Eigen::VectorXd b_eq = Eigen::VectorXd::Zero(2 * (m - r));

        for (size_t i = r; i < m; ++i) {
            const size_t row = 2 * (i - r);
            A_eq(row, i) = 1.0;
            A_eq(row + 1, i) = -1.0;
        }

        A_new = linalg::vstack(A_new, A_eq);
        b_new = linalg::vstack(b_new, b_eq);
    }

    // Eliminate all auxiliary dimensions via Fourier-Motzkin.
    for (size_t i = 0; i < n_aux; ++i) {
        auto eliminated = linalg::fourier_motzkin_elimination(A_new, b_new, m);
        A_new = std::move(eliminated.first);
        b_new = std::move(eliminated.second);
    }

    // Transform constraints from y' to y with y' = U^T y.
    this->A = A_new * U.transpose();
    this->b = std::move(b_new);
}

// --- Binary operations ---

bool do_intersect_impl(const HPolytope& s1, const HPolytope& s2) {
    return !intersection_impl(s1, s2).is_empty();
}

HPolytope intersection_impl(const HPolytope& s1, const HPolytope& s2) {
    return HPolytope(linalg::vstack(s1.A, s2.A), linalg::vstack(s1.b, s2.b), false);
}

HPolytope minkowski_sum_impl(const HPolytope& s1, const HPolytope& s2) {
    const size_t n = s1.dim();
    const size_t h = s1.n_constraints();

    Eigen::MatrixXd A_new = linalg::vstack(
        s1.A, linalg::vstack(Eigen::MatrixXd::Identity(n, n), -Eigen::MatrixXd::Identity(n, n)));

    Eigen::VectorXd b_new = Eigen::VectorXd::Zero(h + 2 * n);
    b_new.head(h) = s1.b;

    // First h constraints: offset with support function of s2.
    for (size_t i = 0; i < h; ++i) {
        Eigen::VectorXd direction = A_new.row(i).transpose();
        b_new(i) += s2.support_function(direction).value;
    }

    // Remaining 2n box constraints: sum supports from both sets.
    for (size_t i = 0; i < 2 * n; ++i) {
        const size_t row = h + i;
        Eigen::VectorXd direction = A_new.row(row).transpose();
        b_new(row) = s1.support_function(direction).value + s2.support_function(direction).value;
    }

    // Merge duplicate constraint directions and keep the tightest bound.
    std::vector<Eigen::Index> keep_rows;
    keep_rows.reserve(static_cast<size_t>(A_new.rows()));

    for (Eigen::Index row = 0; row < A_new.rows(); ++row) {
        bool merged = false;
        for (Eigen::Index kept : keep_rows) {
            if ((A_new.row(row) - A_new.row(kept)).norm() <= 1e-12) {
                if (b_new(row) < b_new(kept)) {
                    b_new(kept) = b_new(row);
                }
                merged = true;
                break;
            }
        }
        if (!merged) {
            keep_rows.push_back(row);
        }
    }

    Eigen::MatrixXd A_reduced(keep_rows.size(), A_new.cols());
    Eigen::VectorXd b_reduced(keep_rows.size());
    for (size_t i = 0; i < keep_rows.size(); ++i) {
        A_reduced.row(i) = A_new.row(keep_rows[i]);
        b_reduced(i) = b_new(keep_rows[i]);
    }

    return HPolytope(std::move(A_reduced), std::move(b_reduced), false);
}

} // namespace geosets
