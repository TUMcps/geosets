#pragma once

#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <cstddef>
#include <cstdio>
#include <geosets/sets/interface.h>
#include <stdexcept>
#include <vector>

namespace geosets {

struct VPolytopeValidationException : std::runtime_error {
    using std::runtime_error::runtime_error;
};

struct VPolytope : public GeometricSet<VPolytope> {
    Eigen::MatrixXd vertices;

  private:
    inline void validate() override {
        if (this->dim() >= this->n_vertices()) {
            throw VPolytopeValidationException(
                "A VPolytope must have more vertices (" + std::to_string(this->n_vertices()) +
                ") than its dimension (" + std::to_string(this->dim()) + ").");
        }
        // TODO: Compute convex hull here!
    }

    // --- GeometricSet implementation hooks ---
    inline bool is_empty_impl() const override { return n_vertices() == 0; }
    bool bounded_impl() const override;
    bool degenerate_impl() const override;
    bool contains_point_impl(const Eigen::VectorXd& point, double tol = 1e-9) const override;
    double volume_impl() const override;
    Eigen::MatrixXd to_vertices_impl() const override;
    Eigen::VectorXd center_impl() const override;
    geosets::Support support_function_impl(const Eigen::VectorXd& direction) const override;
    geosets::Hyperplane tangent_hyperplane_impl(const Eigen::VectorXd& direction) const override;
    Eigen::VectorXd boundary_point_impl(const Eigen::VectorXd& direction) const override;
    std::pair<double, double>
    line_intersection_bounds_impl(const Eigen::VectorXd& point,
                                  const Eigen::VectorXd& direction) const override;
    void translate_impl_(const Eigen::VectorXd& vector) override;
    void linear_transform_impl_(const Eigen::MatrixXd& matrix) override;

  public:
    // --- Constructor ---
    VPolytope() = delete;
    // TODO: implement convex hull
    explicit VPolytope(const Eigen::MatrixXd& vertices, bool validate = true) : vertices(vertices) {
        if (validate) this->validate();
    };
    explicit VPolytope(Eigen::MatrixXd&& vertices, bool validate = true)
        : vertices(std::move(vertices)) {
        if (validate) this->validate();
    };

    explicit VPolytope(const std::vector<Eigen::VectorXd>& vertex_list, bool validate = true) {
        if (vertex_list.empty()) {
            throw VPolytopeValidationException("Vertex list cannot be empty.");
        }
        this->vertices.resize(static_cast<Eigen::Index>(vertex_list.size()), vertex_list[0].size());
        for (size_t i = 0; i < vertex_list.size(); ++i) {
            this->vertices.row(static_cast<Eigen::Index>(i)) = vertex_list[i].transpose();
        }
        if (validate) this->validate();
    }

    static VPolytope make_empty(std::size_t dim) {
        return VPolytope(Eigen::MatrixXd::Zero(0, dim), false);
    }
    static VPolytope from_unit_box(size_t dim);
    static VPolytope from_random(size_t dim, size_t n_vertices);

    void append_vertices(const Eigen::MatrixXd& new_vertices) {
        if (new_vertices.cols() != static_cast<Eigen::Index>(this->dim())) {
            throw VPolytopeValidationException(
                "New vertices must have the same dimension as existing vertices.");
        }
        Eigen::MatrixXd combined(this->vertices.rows() + new_vertices.rows(), this->dim());
        combined << this->vertices, new_vertices;
        this->vertices = std::move(combined);
    }

    std::string str() const override;

    // --- Methods custom to HPolytope ---
    inline size_t n_vertices() const { return vertices.rows(); }

    // --- Simple inlines ---
    inline std::size_t dim() const override { return vertices.cols(); };
};

// --- Binary operations ---
bool do_intersect_impl(const VPolytope& s1, const VPolytope& s2);
VPolytope intersection_impl(const VPolytope& s1, const VPolytope& s2);
VPolytope minkowski_sum_impl(const VPolytope& s1, const VPolytope& s2);

} // namespace geosets
