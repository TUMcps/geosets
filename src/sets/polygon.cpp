#include "geosets/exceptions.h"
#include "geosets/structs.h"
#include <algorithm>
#include <cmath>
#include <geosets/core.h>
#include <geosets/sets/polygon.h>
#include <limits>

namespace geosets {

namespace {

double cross_2d(const Eigen::Vector2d& lhs, const Eigen::Vector2d& rhs) {
    return lhs.x() * rhs.y() - lhs.y() * rhs.x();
}

void update_line_bounds(double t, double atol, double& lower, double& upper) {
    if (t <= atol) lower = std::max(lower, t);
    if (t >= -atol) upper = std::min(upper, t);
}

} // namespace

Polygon Polygon::from_unit_box(size_t dim) {
    if (dim != 2) {
        throw PolygonDimensionException();
    }

    BoostRing ring;
    ring.reserve(5);
    ring.emplace_back(1.0, 1.0);
    ring.emplace_back(1.0, -1.0);
    ring.emplace_back(-1.0, -1.0);
    ring.emplace_back(-1.0, 1.0);
    ring.emplace_back(1.0, 1.0); // Close the ring
    return Polygon(std::move(ring));
}

Polygon Polygon::from_random(size_t dim, size_t n_vertices) {
    if (dim != 2) {
        throw PolygonDimensionException();
    }

    auto points = Eigen::MatrixXd::Random(n_vertices, 2).eval();
    BoostRing ring;
    ring.reserve(static_cast<size_t>(n_vertices + 1));
    for (Eigen::Index i = 0; i < points.rows(); ++i) {
        ring.emplace_back(points(i, 0), points(i, 1));
    }
    ring.emplace_back(points(0, 0), points(0, 1)); // Close the ring

    // Compute centroid
    double cx = 0.0, cy = 0.0;
    for (size_t i = 0; i < ring.size() - 1; ++i) {
        cx += ring[i].x();
        cy += ring[i].y();
    }
    cx /= static_cast<double>(ring.size() - 1);
    cy /= static_cast<double>(ring.size() - 1);

    // Sort points by angle from centroid (clockwise)
    std::sort(ring.begin(), ring.end() - 1, [cx, cy](const BoostPoint& a, const BoostPoint& b) {
        double angle_a = std::atan2(a.y() - cy, a.x() - cx);
        double angle_b = std::atan2(b.y() - cy, b.x() - cx);
        return angle_a > angle_b; // Clockwise ordering
    });

    // Update the closing point
    ring.back() = ring.front();

    return Polygon(std::move(ring));
}

std::string Polygon::str() const {
    std::stringstream ss;
    Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n", "[", "]");
    ss << "Vertices:\n";
    for (const auto& point : this->boost_polygon.outer()) {
        ss << "[" << point.x() << ", " << point.y() << "]\n";
    }
    return "Polygon\n" + ss.str();
}

bool Polygon::degenerate_impl() const {
    return std::abs(boost::geometry::area(boost_polygon)) <= this->atol_;
}

bool Polygon::bounded_impl() const {
    for (const auto& point : this->boost_polygon.outer()) {
        if (!std::isfinite(point.x()) || !std::isfinite(point.y())) {
            return false;
        }
    }
    return true;
}

bool Polygon::contains_point_impl(const Eigen::VectorXd& point, double tol) const {
    BoostPoint p(point(0), point(1));

    if (boost::geometry::covered_by(p, boost_polygon)) {
        return true;
    }

    if (tol <= 0.0) {
        return false;
    }

    return boost::geometry::distance(p, boost_polygon) <= tol;
}

Eigen::MatrixXd Polygon::to_vertices_impl() const {
    Eigen::MatrixXd vertices(this->n_vertices(), 2);

    for (Eigen::Index i = 0; i < vertices.rows(); ++i) {
        const auto& point = this->boost_polygon.outer()[static_cast<size_t>(i)];
        vertices(i, 0) = point.x();
        vertices(i, 1) = point.y();
    }

    return vertices;
}

Eigen::VectorXd Polygon::center_impl() const {
    BoostPoint centroid{0.0, 0.0};
    boost::geometry::centroid(boost_polygon, centroid);

    Eigen::VectorXd c(2);
    c << centroid.x(), centroid.y();
    return c;
}

geosets::Hyperplane Polygon::tangent_hyperplane_impl(const Eigen::VectorXd&) const {
    throw geosets::NotImplemented("tangent_hyperplane is not implemented for Polygon yet.");
}

geosets::Support Polygon::support_function_impl(const Eigen::VectorXd& direction) const {
    const auto& ring = this->boost_polygon.outer();
    if (ring.size() < 2) {
        throw PolygonValidationException("Polygon has insufficient vertices for support function.");
    }

    const double dx = direction(0);
    const double dy = direction(1);

    double best_x = ring[0].x();
    double best_y = ring[0].y();
    double best_value = best_x * dx + best_y * dy;

    for (size_t i = 1; i + 1 < ring.size(); ++i) {
        const double x = ring[i].x();
        const double y = ring[i].y();
        const double value = x * dx + y * dy;
        if (value > best_value) {
            best_value = value;
            best_x = x;
            best_y = y;
        }
    }

    return {best_value, Eigen::Vector2d(best_x, best_y)};
}

double Polygon::volume_impl() const {
    const double area = std::abs(boost::geometry::area(boost_polygon));
    if (area <= this->atol_) {
        return 0.0;
    }
    return area;
}

Eigen::VectorXd Polygon::boundary_point_impl(const Eigen::VectorXd&) const {
    throw geosets::NotImplemented(
        "Boundary point computation for polygons is not yet implemented.");
}

std::pair<double, double>
Polygon::line_intersection_bounds_impl(const Eigen::VectorXd& point,
                                       const Eigen::VectorXd& direction) const {
    const Eigen::Vector2d p(point(0), point(1));
    const Eigen::Vector2d d(direction(0), direction(1));
    const double direction_norm_sq = d.squaredNorm();

    double lower = -std::numeric_limits<double>::infinity();
    double upper = std::numeric_limits<double>::infinity();

    const auto& ring = this->boost_polygon.outer();
    for (size_t i = 0; i + 1 < ring.size(); ++i) {
        const Eigen::Vector2d edge_start(ring[i].x(), ring[i].y());
        const Eigen::Vector2d edge_end(ring[i + 1].x(), ring[i + 1].y());
        const Eigen::Vector2d edge = edge_end - edge_start;
        const Eigen::Vector2d offset = edge_start - p;
        const double denominator = cross_2d(d, edge);

        if (std::abs(denominator) <= this->atol_) {
            if (std::abs(cross_2d(offset, d)) <= this->atol_) {
                update_line_bounds(offset.dot(d) / direction_norm_sq, this->atol_, lower, upper);
                update_line_bounds((edge_end - p).dot(d) / direction_norm_sq, this->atol_, lower,
                                   upper);
            }
            continue;
        }

        const double t = cross_2d(offset, edge) / denominator;
        const double u = cross_2d(offset, d) / denominator;
        if (u >= -this->atol_ && u <= 1.0 + this->atol_) {
            update_line_bounds(t, this->atol_, lower, upper);
        }
    }

    if (!std::isfinite(lower) || !std::isfinite(upper)) {
        throw std::runtime_error("The line does not intersect the Polygon boundary.");
    }
    if (lower > upper + this->atol_) {
        throw std::runtime_error("The line does not intersect the Polygon.");
    }

    return {lower, upper};
}

void Polygon::translate_impl_(const Eigen::VectorXd& vector) {
    const double dx = vector(0);
    const double dy = vector(1);

    for (auto& point : boost_polygon.outer()) {
        point.x(point.x() + dx);
        point.y(point.y() + dy);
    }
}

void Polygon::linear_transform_impl_(const Eigen::MatrixXd& matrix) {
    if (matrix.rows() != 2) {
        throw PolygonDimensionException();
    }

    const double a00 = matrix(0, 0);
    const double a01 = matrix(0, 1);
    const double a10 = matrix(1, 0);
    const double a11 = matrix(1, 1);

    for (auto& point : this->boost_polygon.outer()) {
        const double x = point.x();
        const double y = point.y();
        point.x(a00 * x + a01 * y);
        point.y(a10 * x + a11 * y);
    }

    boost::geometry::correct(this->boost_polygon);
}

Polygon intersection_impl(const Polygon& s1, const Polygon& s2) {
    std::vector<BoostPolygon> output;
    boost::geometry::intersection(s1.boost_polygon, s2.boost_polygon, output);

    if (output.empty()) {
        return Polygon::make_empty(s1.dim());
    }

    const BoostPolygon& result = output.front();
    if (result.outer().size() < 4) {
        return Polygon::make_empty(s1.dim());
    }

    BoostRing ring = result.outer();
    return Polygon(std::move(ring), false);
}

bool do_intersect_impl(const Polygon& s1, const Polygon& s2) {
    return boost::geometry::intersects(s1.boost_polygon, s2.boost_polygon);
}

Polygon minkowski_sum_impl(const Polygon& s1, const Polygon& s2) {
    const auto v1 = s1.to_vertices();
    const auto v2 = s2.to_vertices();

    std::vector<BoostPoint> summed;
    summed.reserve(static_cast<size_t>(v1.rows() * v2.rows()));

    for (Eigen::Index i = 0; i < v1.rows(); ++i) {
        for (Eigen::Index j = 0; j < v2.rows(); ++j) {
            summed.emplace_back(v1(i, 0) + v2(j, 0), v1(i, 1) + v2(j, 1));
        }
    }

    boost::geometry::model::multi_point<BoostPoint> points;
    points.reserve(summed.size());
    for (const auto& point : summed) {
        points.push_back(point);
    }

    BoostPolygon hull;
    boost::geometry::convex_hull(points, hull);

    if (hull.outer().size() < 4) {
        return Polygon::make_empty(s1.dim());
    }

    BoostRing ring = hull.outer();
    return Polygon(std::move(ring), false);
}

Polygon set_difference_impl(const Polygon& s1, const Polygon& s2) {
    std::vector<BoostPolygon> output;
    boost::geometry::difference(s1.boost_polygon, s2.boost_polygon, output);

    if (output.empty()) {
        return Polygon::make_empty(s1.dim());
    }

    return Polygon(std::move(output.front()), false);
}

} // namespace geosets
