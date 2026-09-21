#include "geosets/structs.h"
#include <Eigen/Dense>
#include <cmath>
#include <cstddef>
#include <geosets/conversion/to_hpolytope.h>
#include <geosets/exceptions.h>
#include <geosets/highs_wrapper.h>
#include <geosets/linalg_utils.h>
#include <geosets/qhull_wrapper.h>
#include <geosets/sets/zonotope.h>
#include <numeric>
#include <vector>

namespace geosets {

// --- Static constructors ---

Zonotope Zonotope::from_unit_box(size_t dim) {
    auto c = Eigen::VectorXd::Zero(dim);
    auto G = Eigen::MatrixXd::Identity(dim, dim);
    return Zonotope(std::move(c), std::move(G));
}

Zonotope Zonotope::from_random(size_t dim, size_t n_generators) {
    auto c = Eigen::VectorXd::Random(dim);
    auto G = Eigen::MatrixXd::Random(dim, n_generators);
    return Zonotope(std::move(c), std::move(G));
}

// --- Properties ---

std::string Zonotope::str() const {
    std::stringstream ss;
    Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n", "[", "]");
    ss << "c:\n" << this->c.format(fmt) << "\n";
    ss << "G:\n" << this->G.format(fmt);
    return "Zonotope\n" + ss.str();
}

bool Zonotope::bounded_impl() const {
    return true;
}

bool Zonotope::degenerate_impl() const {
    return linalg::rank(this->G, this->atol_) < this->dim();
}

// --- Core operations ---

Eigen::VectorXd Zonotope::center_impl() const {
    return this->c;
}

geosets::Hyperplane Zonotope::tangent_hyperplane_impl(const Eigen::VectorXd&) const {
    throw geosets::NotImplemented("tangent_hyperplane is not implemented for Zonotope yet.");
}

geosets::Support Zonotope::support_function_impl(const Eigen::VectorXd& direction) const {
    // val = d^T * c + sum(|d^T * g_i|)
    // vec = c + G * sign(G^T * d)
    Eigen::VectorXd Gtd = this->G.transpose() * direction;
    double value = direction.dot(this->c) + Gtd.array().abs().sum();
    Eigen::VectorXd signs = Gtd.unaryExpr([](double x) { return x >= 0.0 ? 1.0 : -1.0; });
    Eigen::VectorXd vector = this->c + this->G * signs;
    return {value, vector};
}

double Zonotope::zonotope_norm(const Eigen::VectorXd& point) const {
    // Minimize t such that G * beta = point, |beta_i| <= t
    // Variables: [t, beta_1, ..., beta_p]  (p+1 variables)
    // Constraints:
    //   G * beta  <= point   (equality via two inequalities)
    //  -G * beta  <= -point
    //   beta_i - t <= 0      (i.e., beta_i <= t)
    //  -beta_i - t <= 0      (i.e., -beta_i <= t)
    const size_t n = this->dim();
    const size_t p = this->n_generators();

    if (p == 0) {
        return point.norm() < this->atol_ ? 0.0 : std::numeric_limits<double>::infinity();
    }

    const size_t n_vars = p + 1; // [t, beta_1, ..., beta_p]
    const size_t m = 2 * n + 2 * p;
    Eigen::MatrixXd A_lp = Eigen::MatrixXd::Zero(m, n_vars);
    Eigen::VectorXd b_lp(m);

    // G * beta <= point
    A_lp.block(0, 1, n, p) = this->G;
    b_lp.head(n) = point;

    // -G * beta <= -point
    A_lp.block(n, 1, n, p) = -this->G;
    b_lp.segment(n, n) = -point;

    // beta_i - t <= 0
    for (size_t i = 0; i < p; ++i) {
        A_lp(2 * n + i, 0) = -1.0;    // -t
        A_lp(2 * n + i, 1 + i) = 1.0; // +beta_i
    }
    b_lp.segment(2 * n, p) = Eigen::VectorXd::Zero(p);

    // -beta_i - t <= 0
    for (size_t i = 0; i < p; ++i) {
        A_lp(2 * n + p + i, 0) = -1.0;     // -t
        A_lp(2 * n + p + i, 1 + i) = -1.0; // -beta_i
    }
    b_lp.tail(p) = Eigen::VectorXd::Zero(p);

    // Objective: minimize t (highs_wrapper maximizes, so we pass t as objective)
    // highs_wrapper negates internally to maximize, so passing -t gives minimization
    Eigen::VectorXd objective = Eigen::VectorXd::Zero(n_vars);
    objective(0) = -1.0; // maximize -t = minimize t

    Eigen::VectorXd solution;
    const HighsModelStatus status = highs_wrapper::solve_lp(A_lp, b_lp, objective, &solution);

    if (highs_wrapper::is_infeasible_status(status)) {
        return std::numeric_limits<double>::infinity();
    }
    if (status != HighsModelStatus::kOptimal) {
        throw FailedOptimizationProblem();
    }

    return solution(0);
}

bool Zonotope::contains_point_impl(const Eigen::VectorXd& point, double tol) const {
    Eigen::VectorXd d = point - this->c;
    Zonotope zero_centered(Eigen::VectorXd::Zero(this->dim()), this->G, false);
    return zero_centered.zonotope_norm(d) <= 1.0 + tol;
}

Eigen::VectorXd Zonotope::boundary_point_impl(const Eigen::VectorXd& direction) const {
    if (direction.norm() < this->atol_) return this->c;

    Zonotope zero_centered(Eigen::VectorXd::Zero(this->dim()), this->G, false);
    double norm = zero_centered.zonotope_norm(direction);

    if (std::isinf(norm)) {
        throw DegenerateSetException();
    }
    if (norm < this->atol_) return this->c;

    return this->c + direction / norm;
}

std::pair<double, double>
Zonotope::line_intersection_bounds_impl(const Eigen::VectorXd& point,
                                        const Eigen::VectorXd& direction) const {

    (void)point;     // for static check
    (void)direction; // for static check
    throw NotImplemented();
}

double Zonotope::volume_impl() const {
    const size_t n = this->dim();
    const size_t p = this->n_generators();

    if (p < n) return 0.0;

    // V = 2^n * sum |det(G(:,S))| over all n-combinations of p generators
    double det_sum = 0.0;

    // Generate all C(p,n) combinations using index array
    std::vector<size_t> indices(n);
    std::iota(indices.begin(), indices.end(), 0);

    while (true) {
        // Extract n x n submatrix for this combination
        Eigen::MatrixXd sub(n, n);
        for (size_t j = 0; j < n; ++j) {
            sub.col(j) = this->G.col(indices[j]);
        }
        det_sum += std::abs(sub.determinant());

        // Advance to next combination
        size_t i = n;
        while (i > 0) {
            --i;
            ++indices[i];
            if (indices[i] <= p - n + i) break;
            if (i == 0) goto done;
        }
        // Fill remaining indices
        for (size_t j = i + 1; j < n; ++j) {
            indices[j] = indices[j - 1] + 1;
        }
    }
done:
    return std::pow(2.0, static_cast<double>(n)) * det_sum;
}

Eigen::MatrixXd Zonotope::to_vertices_impl() const {
    const size_t n = this->dim();
    const size_t p = this->n_generators();

    if (p == 0) {
        // Single point: the center
        Eigen::MatrixXd verts(1, n);
        verts.row(0) = this->c.transpose();
        return verts;
    }

    // Enumerate all 2^p sign combinations to get candidate vertices
    const size_t num_candidates = static_cast<size_t>(1) << p;
    Eigen::MatrixXd candidates(num_candidates, n);

    for (size_t i = 0; i < num_candidates; ++i) {
        Eigen::VectorXd v = this->c;
        for (size_t j = 0; j < p; ++j) {
            if ((i >> j) & 1) {
                v += this->G.col(j);
            } else {
                v -= this->G.col(j);
            }
        }
        candidates.row(i) = v.transpose();
    }

    // Take convex hull to get actual vertices
    if (linalg::rank(this->G, this->atol_) < n) {
        // Degenerate case: return unique candidates
        return candidates;
    }

    auto hull = qhull::build_convex_hull(candidates);
    return qhull::vertices_matrix(hull);
}

void Zonotope::translate_impl_(const Eigen::VectorXd& vector) {
    if (is_empty_impl()) return;
    this->c += vector;
}

void Zonotope::linear_transform_impl_(const Eigen::MatrixXd& matrix) {
    if (is_empty_impl()) return;
    this->c = matrix * this->c;
    this->G = matrix * this->G;
}

// --- Binary operations ---

bool do_intersect_impl(const Zonotope& /*s1*/, const Zonotope& /*s2*/) {
    throw NotImplemented("do_intersect is not implemented for Zonotope.");
}

Zonotope intersection_impl(const Zonotope& /*s1*/, const Zonotope& /*s2*/) {
    throw NotImplemented("intersection is not implemented for Zonotope.");
}

Zonotope minkowski_sum_impl(const Zonotope& s1, const Zonotope& s2) {
    return Zonotope(s1.c + s2.c, linalg::hstack(s1.G, s2.G), false);
}

} // namespace geosets
