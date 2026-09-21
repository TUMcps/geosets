#pragma once

#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <cstddef>
#include <cstdio>
#include <geosets/sets/interface.h>

namespace geosets {

struct IntervalValidationException : std::runtime_error {
    using std::runtime_error::runtime_error;
};

struct Interval : public GeometricSet<Interval> {
    Eigen::VectorXd lb;
    Eigen::VectorXd ub;

  private:
    inline void validate() override {
        if (lb.rows() != ub.rows() || lb.cols() != 1 || ub.cols() != 1) {
            throw IntervalValidationException(
                "Interval: lb and ub must be column vectors of same size");
        }
        if ((ub - lb).minCoeff() < 0.0) {
            throw IntervalValidationException("Interval: ub must be greater than or equal to lb");
        }
    }

    // --- GeometricSet implementation hooks ---
    bool is_empty_impl() const override { return (ub - lb).maxCoeff() <= 0.0; }
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
    Interval() = delete;
    // Values are copied (for lvalues)
    explicit Interval(const Eigen::VectorXd& lb, const Eigen::VectorXd& ub, bool validate = true)
        : lb(lb), ub(ub) {
        if (validate) this->validate();
    }
    // Values are moved (for rvalues)
    explicit Interval(Eigen::VectorXd&& lb, Eigen::VectorXd&& ub, bool validate = true)
        : lb(std::move(lb)), ub(std::move(ub)) {
        if (validate) this->validate();
    }

    static Interval make_empty(std::size_t dim) {
        return Interval(Eigen::VectorXd::Zero(dim), Eigen::VectorXd::Zero(dim), false);
    }
    static Interval from_unit_box(size_t dim);
    static Interval from_random(size_t dim);

    std::string str() const override;

    // --- Simple inlines ---
    inline std::size_t dim() const override { return lb.size(); }
};

// --- Binary operations ---
bool do_intersect_impl(const Interval& a, const Interval& b);
Interval intersection_impl(const Interval& a, const Interval& b);
Interval minkowski_sum_impl(const Interval& a, const Interval& b);

} // namespace geosets
