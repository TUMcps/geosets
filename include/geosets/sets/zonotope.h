#pragma once

#include <Eigen/Dense>
#include <cstddef>
#include <geosets/sets/interface.h>

namespace geosets {

struct ZonotopeValidationException : std::runtime_error {
    using std::runtime_error::runtime_error;
};

struct Zonotope : public GeometricSet<Zonotope> {
    Eigen::VectorXd c;
    Eigen::MatrixXd G;

  private:
    inline void validate() override {
        if (c.size() != G.rows()) {
            throw ZonotopeValidationException("Center size (" + std::to_string(c.size()) +
                                              ") must match generator row count (" +
                                              std::to_string(G.rows()) + ").");
        }
    }

    // --- GeometricSet implementation hooks ---
    inline bool is_empty_impl() const override { return c.size() == 0; }
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
    Zonotope() = delete;
    explicit Zonotope(const Eigen::VectorXd& c, const Eigen::MatrixXd& G, bool validate = true)
        : c(c), G(G) {
        if (validate) this->validate();
    }
    explicit Zonotope(Eigen::VectorXd&& c, Eigen::MatrixXd&& G, bool validate = true)
        : c(std::move(c)), G(std::move(G)) {
        if (validate) this->validate();
    }

    static Zonotope make_empty(std::size_t dim) {
        return Zonotope(Eigen::VectorXd(0), Eigen::MatrixXd(0, dim), false);
    }
    static Zonotope from_unit_box(size_t dim);
    static Zonotope from_random(size_t dim, size_t n_generators);

    std::string str() const override;

    // --- Methods custom to Zonotope ---
    inline size_t n_generators() const { return G.cols(); }
    inline double order() const { return static_cast<double>(n_generators()) / dim(); }

    // Computes the zonotope norm of a point w.r.t. this zonotope.
    // The zonotope must be centered at the origin.
    // Returns the minimum t such that point = G * beta with |beta_i| <= t.
    // Returns infinity if the point is not reachable by any scaling.
    double zonotope_norm(const Eigen::VectorXd& point) const;

    // --- Simple inlines ---
    inline std::size_t dim() const override {
        // For empty sets, c.size() == 0 but G has shape (0, dim)
        return is_empty_impl() ? G.cols() : c.size();
    };
};

// --- Binary operations ---
bool do_intersect_impl(const Zonotope& s1, const Zonotope& s2);
Zonotope intersection_impl(const Zonotope& s1, const Zonotope& s2);
Zonotope minkowski_sum_impl(const Zonotope& s1, const Zonotope& s2);

} // namespace geosets
