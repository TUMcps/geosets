#pragma once

#include <Eigen/Dense>
#include <cstddef>
#include <geosets/exceptions.h>
#include <geosets/structs.h>
#include <stdexcept>
#include <string>
#include <utility>

namespace geosets {

/*
 * GeometricSetInterface is the pure virtual base class
 */
struct GeometricSetInterface {
    virtual ~GeometricSetInterface() = default;
    virtual std::size_t dim() const = 0;
    virtual std::string str() const = 0;

    virtual bool is_empty() const = 0;
    virtual bool is_infinite() const = 0;
    virtual bool bounded() const = 0;
    virtual bool degenerate() const = 0;
    virtual bool contains_point(const Eigen::VectorXd& point, double tol = 1e-9,
                                bool validate = true) const = 0;
    virtual double volume() const = 0;

    /* Returns the set in verterx representation. The vertices are always a copy! */
    virtual Eigen::MatrixXd to_vertices() const = 0;

    virtual Eigen::VectorXd center(bool validate = true) const = 0;
    virtual geosets::Support support_function(const Eigen::VectorXd& direction,
                                              bool validate = true) const = 0;
    virtual geosets::Hyperplane tangent_hyperplane(const Eigen::VectorXd& direction,
                                                   bool validate = true) const = 0;
    virtual Eigen::VectorXd boundary_point(const Eigen::VectorXd& direction,
                                           bool validate = true) const = 0;
    virtual std::pair<double, double> line_intersection_bounds(const Eigen::VectorXd& point,
                                                               const Eigen::VectorXd& direction,
                                                               bool validate = true) const = 0;
    virtual std::pair<Eigen::VectorXd, Eigen::VectorXd>
    line_intersection(const Eigen::VectorXd& point, const Eigen::VectorXd& direction,
                      bool validate = true) const = 0;

    // --- In place methods ---
    // In place needs to be defined by the derived class
    virtual void translate_(const Eigen::VectorXd& vector) = 0;
    virtual void linear_transform_(const Eigen::MatrixXd& matrix) = 0;

    // --- Purely virtual pointer-based methods ---
    // These are needed for the python interface to work
    // virtual GeometricSetInterface* translate_ptr(const Eigen::VectorXd& v) const = 0;
    // virtual GeometricSetInterface* linear_transform_ptr(const Eigen::MatrixXd& m) const = 0;
};

/*
 * GeometricSet is the CRTP base class.
 */
template <typename Derived> class GeometricSet : public GeometricSetInterface {
  public:
    virtual ~GeometricSet() = default;

    inline Derived copy() const { return Derived(*static_cast<const Derived*>(this)); }
    inline Derived* copy_to_ptr() const { return new Derived(*static_cast<const Derived*>(this)); }

    inline bool is_infinite() const override { return this->infinite_; }

    static Derived make_infinite(std::size_t dim) {
        Derived obj = Derived::make_empty(dim);
        obj.infinite_ = true;
        return obj;
    }

    // --- Public interface methods with infinite and other checks ---
    // These all call the _impl versions of the derived classes

    bool is_empty() const override {
        if (infinite_) return false;
        return this->is_empty_impl();
    }

    bool bounded() const override {
        if (infinite_) return false;
        if (this->is_empty()) return true;
        return this->bounded_impl();
    }

    bool degenerate() const override {
        if (infinite_) return false;
        return this->degenerate_impl();
    }

    bool contains_point(const Eigen::VectorXd& point, double tol = 1e-9,
                        bool validate = true) const override {
        const double tol_abs = std::abs(tol);
        if (validate) {
            if (infinite_) return true;
            if (this->is_empty()) return false; // potentially redundant
            check_operand_dim(point.size());
        }
        return this->contains_point_impl(point, tol_abs);
    }

    double volume() const override {
        if (infinite_) return std::numeric_limits<double>::infinity();
        if (this->is_empty() || this->degenerate()) return 0.0;
        return this->volume_impl();
    }

    Eigen::MatrixXd to_vertices() const override {
        if (infinite_) throw InfiniteSetException();
        if (this->is_empty()) throw EmptySetException();
        return this->to_vertices_impl();
    }

    Eigen::VectorXd center(bool validate = true) const override {
        if (validate) {
            if (infinite_) throw InfiniteSetException();
            if (this->is_empty()) throw EmptySetException();
        }
        return this->center_impl();
    }

    geosets::Support support_function(const Eigen::VectorXd& direction,
                                      bool validate = true) const override {
        if (validate) {
            if (infinite_) {
                // Note: technically the support vector is undefined here.
                geosets::Support sup;
                sup.value = std::numeric_limits<double>::infinity();
                sup.vector = Eigen::VectorXd::Zero(this->dim());
                return sup;
            }
            if (this->is_empty()) {
                geosets::Support sup;
                sup.value = -std::numeric_limits<double>::infinity();
                sup.vector = Eigen::VectorXd::Zero(this->dim());
                return sup;
            }
            check_operand_dim(direction.size());
        }
        return this->support_function_impl(direction);
    }

    geosets::Hyperplane tangent_hyperplane(const Eigen::VectorXd& direction,
                                           bool validate = true) const override {
        if (validate) {
            // if (infinite_) {
            //     geosets::Hyperplane hp;
            //     hp.b = std::numeric_limits<double>::infinity();
            //     hp.a = Eigen::VectorXd::Zero(this->dim());
            //     return hp;
            // }
            // if (this->is_empty()) {
            //     geosets::Hyperplane hp;
            //     hp.b = -std::numeric_limits<double>::infinity();
            //     hp.a = Eigen::VectorXd::Zero(this->dim());
            //     return hp;
            // }
            if (infinite_) throw InfiniteSetException();
            if (this->is_empty()) throw EmptySetException();
            check_operand_dim(direction.size());
        }
        return this->tangent_hyperplane_impl(direction);
    }

    Eigen::VectorXd boundary_point(const Eigen::VectorXd& direction,
                                   bool validate = true) const override {
        if (validate) {
            if (infinite_) throw InfiniteSetException();
            if (this->is_empty()) throw EmptySetException();
            check_operand_dim(direction.size());
        }
        return this->boundary_point_impl(direction);
    }

    std::pair<double, double> line_intersection_bounds(const Eigen::VectorXd& point,
                                                       const Eigen::VectorXd& direction,
                                                       bool validate = true) const override {
        if (validate) {
            check_operand_dim(point.size());
            check_operand_dim(direction.size());
            if (!point.array().isFinite().all() || !direction.array().isFinite().all()) {
                throw std::invalid_argument("Line point and direction must contain finite values.");
            }
            if (direction.norm() <= this->atol_) {
                throw std::invalid_argument("Line direction must be non-zero.");
            }
            if (infinite_) throw InfiniteSetException();
            if (this->is_empty()) throw EmptySetException();
            if (!this->contains_point_impl(point, this->atol_)) {
                throw std::invalid_argument("Line point must be contained in the set.");
            }
        }
        return this->line_intersection_bounds_impl(point, direction);
    }

    std::pair<Eigen::VectorXd, Eigen::VectorXd>
    line_intersection(const Eigen::VectorXd& point, const Eigen::VectorXd& direction,
                      bool validate = true) const override {
        const auto [lower, upper] = this->line_intersection_bounds(point, direction, validate);
        return {point + lower * direction, point + upper * direction};
    }

    void translate_(const Eigen::VectorXd& vector) override {
        if (infinite_) return;
        check_operand_dim(vector.size());
        this->translate_impl_(vector);
    }

    void linear_transform_(const Eigen::MatrixXd& matrix) override {
        if (infinite_) return;
        this->check_operand_dim(matrix.cols());
        this->linear_transform_impl_(matrix);
    }

    // -- Copy methods just call the in place methods --
    Derived translate(const Eigen::VectorXd& vector) const {
        Derived copy = this->copy();
        copy.translate_(vector);
        return copy;
    }
    Derived linear_transform(const Eigen::MatrixXd& matrix) const {
        Derived copy = this->copy();
        copy.linear_transform_(matrix);
        return copy;
    }

    // Pointer-based methods implementation
    Derived* translate_ptr(const Eigen::VectorXd& v) const {
        return new Derived(this->translate(v));
    }
    Derived* linear_transform_ptr(const Eigen::MatrixXd& m) const {
        return new Derived(this->linear_transform(m));
    }

  protected:
    // Commonly used absolute tolerance for operations
    double atol_ = 1e-9;
    bool infinite_ = false;

    inline void check_operand_dim(size_t other_dim) const {
        size_t this_dim = this->dim();
        if (this_dim != other_dim) {
            throw DimensionMismatchException(this_dim, other_dim);
        }
    }
    virtual void validate() = 0;

    // --- Pure virtual implementation methods ---
    virtual bool is_empty_impl() const = 0;
    virtual bool bounded_impl() const = 0;
    virtual bool degenerate_impl() const = 0;
    virtual bool contains_point_impl(const Eigen::VectorXd& point, double tol) const = 0;
    virtual double volume_impl() const = 0;
    virtual Eigen::MatrixXd to_vertices_impl() const = 0;
    virtual Eigen::VectorXd center_impl() const = 0;
    virtual geosets::Support support_function_impl(const Eigen::VectorXd& direction) const = 0;
    virtual geosets::Hyperplane tangent_hyperplane_impl(const Eigen::VectorXd& direction) const = 0;
    virtual Eigen::VectorXd boundary_point_impl(const Eigen::VectorXd& direction) const = 0;
    virtual std::pair<double, double>
    line_intersection_bounds_impl(const Eigen::VectorXd&, const Eigen::VectorXd&) const = 0;
    virtual void translate_impl_(const Eigen::VectorXd& vector) = 0;
    virtual void linear_transform_impl_(const Eigen::MatrixXd& matrix) = 0;
};

} // namespace geosets
