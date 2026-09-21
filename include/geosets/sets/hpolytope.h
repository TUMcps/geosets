#pragma once

#include "geosets/structs.h"
#include <Eigen/Dense>
#include <cstddef>
#include <cstdio>
#include <geosets/sets/interface.h>
#include <optional>

namespace geosets {

struct HPolytopeValidationException : std::runtime_error {
    using std::runtime_error::runtime_error;
};

// Note: we should probably use libcdd for some operations like redundancy reduction, conversion.

struct HPolytope : public GeometricSet<HPolytope> {
    Eigen::MatrixXd A;
    Eigen::VectorXd b;

  private:
    mutable std::optional<Eigen::VectorXd> cached_center_;

    inline void validate() override {
        if (A.rows() != b.size()) {
            throw HPolytopeValidationException("Number of rows of A (" + std::to_string(A.rows()) +
                                               ") must match size of b (" +
                                               std::to_string(b.size()) + ").");
        }
    }

    // --- GeometricSet implementation hooks ---
    bool degenerate_impl() const override;
    bool is_empty_impl() const override;
    bool bounded_impl() const override;
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
    HPolytope() = delete;
    // Values are copied (for lvalues)
    explicit HPolytope(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, bool validate = true)
        : A(A), b(b) {
        if (validate) this->validate();
    }
    // Values are moved (for rvalues)
    explicit HPolytope(Eigen::MatrixXd&& A, Eigen::VectorXd&& b, bool validate = true)
        : A(std::move(A)), b(std::move(b)) {
        if (validate) this->validate();
    }

    static HPolytope make_empty(std::size_t dim) {
        // A HPolytope without constraints is not empty, but R^d!
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(1, dim);
        Eigen::VectorXd b = Eigen::VectorXd::Constant(1, -1.0);
        return HPolytope(A, b, false);
    }
    static HPolytope from_unit_box(size_t dim);
    static HPolytope from_random(size_t dim, size_t n_constraints);

    std::string str() const override;

    // --- Methods custom to HPolytope ---
    inline size_t n_constraints() const { return A.rows(); }

    // --- Simple inlines ---
    inline std::size_t dim() const override { return A.cols(); };
};

// --- Binary operations ---
bool do_intersect_impl(const HPolytope& s1, const HPolytope& s2);
HPolytope intersection_impl(const HPolytope& s1, const HPolytope& s2);
HPolytope minkowski_sum_impl(const HPolytope& s1, const HPolytope& s2);

} // namespace geosets
