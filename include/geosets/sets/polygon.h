#pragma once

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <geosets/sets/interface.h>
#include <stdexcept>

namespace geosets {

using BoostPoint = boost::geometry::model::d2::point_xy<double>;
using BoostRing = boost::geometry::model::ring<BoostPoint>;
using BoostPolygon = boost::geometry::model::polygon<BoostPoint>;

struct PolygonValidationException : std::runtime_error {
    using std::runtime_error::runtime_error;
};

struct PolygonDimensionException : std::runtime_error {
    PolygonDimensionException() : std::runtime_error("Polygon is only defined for 2D.") {}
};

/**
 * @brief Polygon set implementation using Boost.Geometry.
 * The representation is lossy, as it only represents connected Polygons without holes.
 *
 */
struct Polygon : public GeometricSet<Polygon> {
    BoostPolygon boost_polygon;

  private:
    inline void validate() override {
        // TODO: Move to validate?
        boost::geometry::correct(boost_polygon);

        if (!boost::geometry::is_valid(boost_polygon)) {
            throw PolygonValidationException("Invalid polygon geometry!");
        }
    }

    // --- GeometricSet implementation hooks ---
    inline bool is_empty_impl() const override {
        return this->boost_polygon.outer().size() == 0 ||
               boost::geometry::is_empty(this->boost_polygon);
    };
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
    Polygon() = delete;
    explicit Polygon(const Eigen::MatrixXd& vertices, bool validate = true) {
        auto& ring = boost_polygon.outer();
        ring.reserve(static_cast<size_t>(vertices.rows()));
        for (Eigen::Index i = 0; i < vertices.rows(); ++i) {
            ring.emplace_back(vertices(i, 0), vertices(i, 1));
        }
        if (validate) this->validate();
    }

    explicit Polygon(BoostRing&& points, bool validate = true) {
        boost_polygon.outer() = std::move(points);

        if (validate) this->validate();
    }

    explicit Polygon(BoostPolygon&& polygon, bool validate = true)
        : boost_polygon(std::move(polygon)) {
        if (validate) this->validate();
    }

    static Polygon make_empty(std::size_t dim) {
        if (dim != 2) {
            throw PolygonDimensionException();
        }
        return Polygon(Eigen::MatrixXd::Zero(0, 2), false);
    }
    static Polygon from_unit_box(size_t dim);
    static Polygon from_random(size_t dim, size_t n_vertices);

    std::string str() const override;

    // --- Methods custom to HPolytope ---
    inline size_t n_vertices() const { return this->boost_polygon.outer().size() - 1; }

    // --- Simple inlines ---
    inline std::size_t dim() const override { return 2; };
};

// --- Binary operations ---
bool do_intersect_impl(const Polygon& s1, const Polygon& s2);
Polygon intersection_impl(const Polygon& s1, const Polygon& s2);
Polygon minkowski_sum_impl(const Polygon& s1, const Polygon& s2);
Polygon set_difference_impl(const Polygon& s1, const Polygon& s2);

} // namespace geosets
