#include "geosets/conversion/to_hpolytope.h"
#include "geosets/conversion/to_vpolytope.h"
#include <Highs.h>
#include <algorithm>
#include <cmath>
#include <geosets/core.h>
#include <geosets/linalg_utils.h>
#include <geosets/qhull_wrapper.h>
#include <geosets/sets/vpolytope.h>
#include <vector>

namespace geosets {

namespace {

Eigen::MatrixXd ordered_2d_hull_vertices(const Eigen::MatrixXd& vertices, double atol) {
    auto hull = qhull::build_convex_hull(vertices);
    Eigen::MatrixXd hull_vertices = qhull::vertices_matrix(hull);

    const Eigen::Vector2d centroid = hull_vertices.colwise().mean();
    std::vector<int> order(static_cast<size_t>(hull_vertices.rows()));
    for (int i = 0; i < hull_vertices.rows(); ++i) {
        order[static_cast<size_t>(i)] = i;
    }

    std::sort(order.begin(), order.end(), [&](const int lhs, const int rhs) {
        const double angle_l =
            std::atan2(hull_vertices(lhs, 1) - centroid(1), hull_vertices(lhs, 0) - centroid(0));
        const double angle_r =
            std::atan2(hull_vertices(rhs, 1) - centroid(1), hull_vertices(rhs, 0) - centroid(0));
        return angle_l < angle_r;
    });

    std::vector<Eigen::Vector2d> unique_points;
    unique_points.reserve(order.size());
    for (const int idx : order) {
        const Eigen::Vector2d point = hull_vertices.row(idx).transpose();
        bool duplicate = false;
        for (const auto& existing : unique_points) {
            if ((existing - point).norm() <= atol) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) {
            unique_points.push_back(point);
        }
    }

    Eigen::MatrixXd ordered(static_cast<Eigen::Index>(unique_points.size()), 2);
    for (Eigen::Index i = 0; i < ordered.rows(); ++i) {
        ordered.row(i) = unique_points[static_cast<size_t>(i)].transpose();
    }

    return ordered;
}

double cross_2d(const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
    return a.x() * b.y() - a.y() * b.x();
}

bool is_inside_halfplane(const Eigen::Vector2d& p, const Eigen::Vector2d& edge_start,
                         const Eigen::Vector2d& edge_end, double atol) {
    const Eigen::Vector2d edge = edge_end - edge_start;
    const Eigen::Vector2d rel = p - edge_start;
    return cross_2d(edge, rel) >= -atol;
}

Eigen::Vector2d line_intersection(const Eigen::Vector2d& s, const Eigen::Vector2d& e,
                                  const Eigen::Vector2d& a, const Eigen::Vector2d& b, double atol) {
    const Eigen::Vector2d se = e - s;
    const Eigen::Vector2d ab = b - a;
    const double denom = cross_2d(se, ab);

    if (std::abs(denom) <= atol) {
        return e;
    }

    const Eigen::Vector2d as = a - s;
    const double t = cross_2d(as, ab) / denom;
    return s + t * se;
}

} // namespace

// --- Static constructors ---

VPolytope VPolytope::from_unit_box(size_t dim) {
    // Generate all 2^dim corners of the unit box [-1, 1]^dim
    const size_t n_vertices = 1 << dim; // 2^dim
    Eigen::MatrixXd vertices = Eigen::MatrixXd::Zero(n_vertices, dim);

    for (size_t i = 0; i < n_vertices; ++i) {
        for (size_t d = 0; d < dim; ++d) {
            vertices(i, d) = (i & (1 << d)) ? 1.0 : -1.0;
        }
    }
    return VPolytope(std::move(vertices));
}

VPolytope VPolytope::from_random(size_t dim, size_t n_vertices) {
    auto hull = qhull::build_convex_hull(Eigen::MatrixXd::Random(n_vertices, dim).eval());
    auto vertices = qhull::vertices_matrix(hull);
    return VPolytope(std::move(vertices));
}

// --- Properties ---

std::string VPolytope::str() const {
    std::stringstream ss;
    Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n", "[", "]");
    ss << "Vertices:\n" << this->vertices.format(fmt);
    return "VPolytope\n" + ss.str();
}

bool VPolytope::degenerate_impl() const {
    // This can be checked by computing the rank of the centered vertices matrix
    if (this->n_vertices() <= this->dim()) {
        return true;
    }
    Eigen::VectorXd mean = this->vertices.colwise().mean();
    Eigen::MatrixXd centered = this->vertices.rowwise() - mean.transpose();
    return geosets::linalg::rank(centered, this->atol_) < this->dim();
}

bool VPolytope::bounded_impl() const {
    return this->vertices.allFinite();
}

// --- Rest ---

bool VPolytope::contains_point_impl(const Eigen::VectorXd& point, double tol) const {
    // We solve: minimize 0^T * lambda
    // subject to: V^T * lambda = point
    //             1^T * lambda = 1
    //             lambda >= 0
    const size_t n = this->n_vertices();
    const size_t dim = this->dim();

    Highs highs;
    highs.setOptionValue("output_flag", false);
    highs.setOptionValue("log_to_console", false);

    HighsModel lp;
    lp.lp_.num_col_ = n;
    lp.lp_.num_row_ = dim + 1;

    // Objective: minimize 0 (feasibility problem)
    lp.lp_.sense_ = ObjSense::kMinimize;
    lp.lp_.offset_ = 0.0;
    lp.lp_.col_cost_.assign(n, 0.0);
    lp.lp_.col_lower_.assign(n, 0.0); // lambda >= 0
    lp.lp_.col_upper_.assign(n, kHighsInf);

    // Constraints: V^T * lambda = point and sum(lambda) = 1
    lp.lp_.row_lower_.resize(dim + 1);
    lp.lp_.row_upper_.resize(dim + 1);

    // First dim constraints: V^T * lambda ~= point
    for (size_t i = 0; i < dim; ++i) {
        lp.lp_.row_lower_[i] = point(i) - tol;
        lp.lp_.row_upper_[i] = point(i) + tol;
    }

    // Last constraint: sum(lambda) ~= 1
    lp.lp_.row_lower_[dim] = 1.0 - tol;
    lp.lp_.row_upper_[dim] = 1.0 + tol;

    // Build constraint matrix in COO format
    lp.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.lp_.a_matrix_.start_.resize(n + 1);

    int nnz = 0;
    for (size_t j = 0; j < n; ++j) {
        lp.lp_.a_matrix_.start_[j] = nnz;
        // Add V^T entries (vertices)
        for (size_t i = 0; i < dim; ++i) {
            lp.lp_.a_matrix_.index_.push_back(i);
            lp.lp_.a_matrix_.value_.push_back(this->vertices(j, i));
            nnz++;
        }
        // Add 1 for sum constraint
        lp.lp_.a_matrix_.index_.push_back(dim);
        lp.lp_.a_matrix_.value_.push_back(1.0);
        nnz++;
    }
    lp.lp_.a_matrix_.start_[n] = nnz;

    HighsStatus status = highs.passModel(lp);
    if (status != HighsStatus::kOk) {
        return false;
    }

    status = highs.run();
    if (status != HighsStatus::kOk) {
        return false;
    }

    const HighsModelStatus& model_status = highs.getModelStatus();
    return model_status == HighsModelStatus::kOptimal;
}

Eigen::MatrixXd VPolytope::to_vertices_impl() const {
    return this->vertices;
}

/*
 * Returns the centroid of the vertices.
 */
Eigen::VectorXd VPolytope::center_impl() const {
    return this->vertices.colwise().mean();
}

geosets::Hyperplane VPolytope::tangent_hyperplane_impl(const Eigen::VectorXd&) const {
    throw geosets::NotImplemented("tangent_hyperplane is not implemented for VPolytope yet.");
}

geosets::Support VPolytope::support_function_impl(const Eigen::VectorXd& direction) const {
    // max{ d^T * v_i } over all vertices v_i
    Eigen::VectorXd dot_products = this->vertices * direction;
    Eigen::Index max_idx;
    double max_value = dot_products.maxCoeff(&max_idx);
    return {max_value, this->vertices.row(max_idx).transpose()};
}

double VPolytope::volume_impl() const {
    auto hull = qhull::build_convex_hull(this->vertices);
    return qhull::volume(hull);
}

Eigen::VectorXd VPolytope::boundary_point_impl(const Eigen::VectorXd& direction) const {
    auto hpoly = conversion::vpolytope_to_hpolytope(*this);
    return hpoly.boundary_point(direction);
}

std::pair<double, double>
VPolytope::line_intersection_bounds_impl(const Eigen::VectorXd& point,
                                         const Eigen::VectorXd& direction) const {
    (void)point;     // for static check
    (void)direction; // for static check
    throw NotImplemented();
}

void VPolytope::translate_impl_(const Eigen::VectorXd& vector) {
    this->vertices.rowwise() += vector.transpose();
}

void VPolytope::linear_transform_impl_(const Eigen::MatrixXd& matrix) {
    this->vertices = this->vertices * matrix.transpose();
}

// --- Binary operations ---

bool do_intersect_impl(const VPolytope& s1, const VPolytope& s2) {
    // TODO: Maybe there is a more efficient algorithm than just converting
    HPolytope h1 = conversion::vpolytope_to_hpolytope(s1);
    HPolytope h2 = conversion::vpolytope_to_hpolytope(s2);
    HPolytope h_intersection = intersection_impl(h1, h2);
    return !h_intersection.is_empty();
}

VPolytope intersection_impl(const VPolytope& s1, const VPolytope& s2) {
    if (s1.dim() >= 2) {
        // Convert to HPolytope and back
        HPolytope h1 = conversion::vpolytope_to_hpolytope(s1);
        HPolytope h2 = conversion::vpolytope_to_hpolytope(s2);
        HPolytope h_intersection = intersection_impl(h1, h2);
        return conversion::hpolytope_to_vpolytope(h_intersection);
    }

    const double atol = 1e-9;

    Eigen::MatrixXd subject = ordered_2d_hull_vertices(s1.vertices, atol);
    Eigen::MatrixXd clip = ordered_2d_hull_vertices(s2.vertices, atol);

    std::vector<Eigen::Vector2d> output;
    output.reserve(static_cast<size_t>(subject.rows()));
    for (Eigen::Index i = 0; i < subject.rows(); ++i) {
        output.push_back(subject.row(i).transpose());
    }

    for (Eigen::Index i = 0; i < clip.rows(); ++i) {
        const Eigen::Vector2d edge_start = clip.row(i).transpose();
        const Eigen::Vector2d edge_end = clip.row((i + 1) % clip.rows()).transpose();

        if (output.empty()) {
            break;
        }

        std::vector<Eigen::Vector2d> input = output;
        output.clear();

        Eigen::Vector2d s = input.back();
        for (const Eigen::Vector2d& e : input) {
            const bool e_inside = is_inside_halfplane(e, edge_start, edge_end, atol);
            const bool s_inside = is_inside_halfplane(s, edge_start, edge_end, atol);

            if (e_inside) {
                if (!s_inside) {
                    output.push_back(line_intersection(s, e, edge_start, edge_end, atol));
                }
                output.push_back(e);
            } else if (s_inside) {
                output.push_back(line_intersection(s, e, edge_start, edge_end, atol));
            }

            s = e;
        }
    }

    std::vector<Eigen::Vector2d> deduplicated;
    deduplicated.reserve(output.size());
    for (const auto& p : output) {
        bool duplicate = false;
        for (const auto& q : deduplicated) {
            if ((p - q).norm() <= atol) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) {
            deduplicated.push_back(p);
        }
    }

    if (deduplicated.size() < 3) {
        return VPolytope::make_empty(s1.dim());
    }

    Eigen::MatrixXd result_vertices(static_cast<Eigen::Index>(deduplicated.size()), 2);
    for (Eigen::Index i = 0; i < result_vertices.rows(); ++i) {
        result_vertices.row(i) = deduplicated[static_cast<size_t>(i)].transpose();
    }

    auto result_hull = qhull::build_convex_hull(result_vertices);
    Eigen::MatrixXd hull_vertices = qhull::vertices_matrix(result_hull);
    return VPolytope(std::move(hull_vertices), false);
}

VPolytope minkowski_sum_impl(const VPolytope& s1, const VPolytope& s2) {
    // Minkowski sum of two V-polytopes is the convex hull of all pairwise sums of vertices
    const size_t n1 = s1.n_vertices();
    const size_t n2 = s2.n_vertices();
    const size_t dim = s1.dim();

    Eigen::MatrixXd sum_vertices(n1 * n2, dim);
    size_t idx = 0;

    for (size_t i = 0; i < n1; ++i) {
        for (size_t j = 0; j < n2; ++j) {
            sum_vertices.row(idx) = s1.vertices.row(i) + s2.vertices.row(j);
            idx++;
        }
    }

    // Compute convex hull of the pairwise sums
    auto hull = qhull::build_convex_hull(sum_vertices);
    Eigen::MatrixXd hull_vertices = qhull::vertices_matrix(hull);
    return VPolytope(std::move(hull_vertices), false);
}

} // namespace geosets
