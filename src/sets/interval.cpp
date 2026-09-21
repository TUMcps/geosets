#include "geosets/structs.h"
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <Highs.h>
#include <algorithm>
#include <cstdio>
#include <geosets/core.h>
#include <geosets/exceptions.h>
#include <geosets/linalg_utils.h>
#include <geosets/sets/interval.h>
#include <limits>

namespace geosets {

// --- Static constructors ---

Interval Interval::from_unit_box(size_t dim) {
    return Interval(-Eigen::VectorXd::Ones(dim), Eigen::VectorXd::Ones(dim));
}

Interval Interval::from_random(size_t dim) {
    auto lb = Eigen::VectorXd::Random(dim) * 0.5 - Eigen::VectorXd::Ones(dim) * 0.5;
    auto ub = Eigen::VectorXd::Random(dim) * 0.5 + Eigen::VectorXd::Ones(dim) * 0.5;
    return Interval(std::move(lb), std::move(ub));
}

// --- Properties ---

std::string Interval::str() const {
    std::stringstream ss;
    Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n", "[", "]");
    ss << "lb:\n" << this->lb.format(fmt) << "\n";
    ss << "ub:\n" << this->ub.format(fmt);
    return "Interval\n" + ss.str();
}

bool Interval::degenerate_impl() const {
    return ((this->ub - this->lb).array().abs() < this->atol_).any();
}

bool Interval::bounded_impl() const {
    return this->lb.allFinite() && this->ub.allFinite();
}

// --- Rest ---

bool Interval::contains_point_impl(const Eigen::VectorXd& point, double tol) const {
    // Check if point is within the bounds: lb <= point <= ub
    return (point.array() >= (lb.array() - tol)).all() &&
           (point.array() <= (ub.array() + tol)).all();
}

double Interval::volume_impl() const {
    return (ub - lb).prod();
}

Eigen::MatrixXd Interval::to_vertices_impl() const {
    const int n = static_cast<int>(dim());
    const int n_vertices = 1 << n;

    Eigen::MatrixXd vertices = Eigen::MatrixXd::Zero(n_vertices, n);

    for (int i = 0; i < n_vertices; ++i) {
        for (int d = 0; d < n; ++d) {
            vertices(i, d) = (i & (1 << d)) ? this->ub(d) : this->lb(d);
        }
    }

    return vertices;
}

Eigen::VectorXd Interval::center_impl() const {
    return 0.5 * (this->lb + this->ub);
}

geosets::Hyperplane Interval::tangent_hyperplane_impl(const Eigen::VectorXd& direction) const {
    const int n = static_cast<int>(dim());
    Eigen::VectorXd c = this->center_impl();

    double min_t = std::numeric_limits<double>::infinity();
    int bound_type = -1; // 0 = lower, 1 = upper
    int min_dim = -1;

    for (int i = 0; i < n; ++i) {
        if (direction(i) < -this->atol_) {
            double t = (lb(i) - c(i)) / direction(i);
            if (t < min_t) {
                min_t = t;
                bound_type = 0;
                min_dim = i;
            }
        }
        if (direction(i) > this->atol_) {
            double t = (ub(i) - c(i)) / direction(i);
            if (t < min_t) {
                min_t = t;
                bound_type = 1;
                min_dim = i;
            }
        }
    }

    Eigen::VectorXd a = Eigen::VectorXd::Zero(n);
    double b;
    if (bound_type == 0) {
        a(min_dim) = -1.0;
        b = -lb(min_dim);
    } else {
        a(min_dim) = 1.0;
        b = ub(min_dim);
    }

    return {std::move(a), b};
}

geosets::Support Interval::support_function_impl(const Eigen::VectorXd& direction) const {
    // For an interval [lb, ub], the support function in direction d is:
    // max{ d^T x | lb <= x <= ub } = d^T x* where x*_i = ub_i if d_i >= 0, lb_i if d_i < 0

    Eigen::VectorXd support_vector = Eigen::VectorXd::Zero(direction.size());
    double support_value = 0.0;

    for (int i = 0; i < direction.size(); ++i) {
        if (direction(i) >= 0) {
            support_vector(i) = ub(i);
            support_value += direction(i) * ub(i);
        } else {
            support_vector(i) = lb(i);
            support_value += direction(i) * lb(i);
        }
    }

    return {support_value, support_vector};
}

Eigen::VectorXd Interval::boundary_point_impl(const Eigen::VectorXd& direction) const {
    Eigen::VectorXd c = this->center_impl();
    double t_min = std::numeric_limits<double>::infinity();

    for (int i = 0; i < direction.size(); ++i) {
        if (direction(i) > 0) {
            t_min = std::min(t_min, (ub(i) - c(i)) / direction(i));
        } else if (direction(i) < 0) {
            t_min = std::min(t_min, (lb(i) - c(i)) / direction(i));
        }
    }

    if (std::isinf(t_min)) return c;
    return c + t_min * direction;
}

std::pair<double, double>
Interval::line_intersection_bounds_impl(const Eigen::VectorXd& point,
                                        const Eigen::VectorXd& direction) const {
    double lower = -std::numeric_limits<double>::infinity();
    double upper = std::numeric_limits<double>::infinity();

    for (Eigen::Index i = 0; i < direction.size(); ++i) {
        const double coefficient = direction(i);
        if (std::abs(coefficient) <= this->atol_) {
            if (point(i) < this->lb(i) - this->atol_ || point(i) > this->ub(i) + this->atol_) {
                throw std::runtime_error("The line does not intersect the Interval.");
            }
            continue;
        }

        double t_lb = (this->lb(i) - point(i)) / coefficient;
        double t_ub = (this->ub(i) - point(i)) / coefficient;
        if (t_lb > t_ub) std::swap(t_lb, t_ub);

        lower = std::max(lower, t_lb);
        upper = std::min(upper, t_ub);
    }

    if (!std::isfinite(lower) || !std::isfinite(upper)) {
        throw std::runtime_error(
            "The line intersection is unbounded and has fewer than two boundary points.");
    }
    if (lower > upper + this->atol_) {
        throw std::runtime_error("The line does not intersect the Interval.");
    }

    return {lower, upper};
}

void Interval::translate_impl_(const Eigen::VectorXd& vector) {
    this->lb += vector;
    this->ub += vector;
}

void Interval::linear_transform_impl_(const Eigen::MatrixXd& matrix) {
    Eigen::VectorXd mid = (ub + lb) * 0.5;
    Eigen::VectorXd rad = (ub - lb) * 0.5;

    Eigen::VectorXd new_mid = matrix * mid;
    Eigen::VectorXd new_rad = matrix.array().abs().matrix() * rad;

    lb = new_mid - new_rad;
    ub = new_mid + new_rad;
}

// --- Binary operations ---

bool do_intersect_impl(const Interval& s1, const Interval& s2) {
    return !intersection_impl(s1, s2).is_empty();
}

Interval intersection_impl(const Interval& s1, const Interval& s2) {
    // For intervals, intersection is element-wise max of lower bounds and min of upper bounds
    Eigen::VectorXd new_lb = s1.lb.array().max(s2.lb.array());
    Eigen::VectorXd new_ub = s1.ub.array().min(s2.ub.array());

    // Check if the intersection is valid (non-empty)
    if ((new_ub - new_lb).minCoeff() < -1e-9) {
        return Interval::make_empty(s1.dim());
    }

    return Interval(std::move(new_lb), std::move(new_ub), false);
}

Interval minkowski_sum_impl(const Interval& s1, const Interval& s2) {
    // For intervals, Minkowski sum is element-wise sum of bounds
    Eigen::VectorXd new_lb = s1.lb + s2.lb;
    Eigen::VectorXd new_ub = s1.ub + s2.ub;

    return Interval(std::move(new_lb), std::move(new_ub), false);
}

} // namespace geosets
