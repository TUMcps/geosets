#include <boost/geometry.hpp>
#include <boost/geometry/geometries/multi_point.hpp>
#include <geosets/conversion/to_polygon.h>

namespace geosets::conversion {

Polygon vpolytope_to_polygon(const VPolytope& vpolytope) {
    if (vpolytope.dim() != 2) {
        throw PolygonDimensionException();
    }
    if (vpolytope.is_empty()) {
        return Polygon::make_empty(2);
    }

    boost::geometry::model::multi_point<BoostPoint> points;
    points.reserve(static_cast<size_t>(vpolytope.n_vertices()));
    for (Eigen::Index i = 0; i < vpolytope.vertices.rows(); ++i) {
        points.push_back(BoostPoint(vpolytope.vertices(i, 0), vpolytope.vertices(i, 1)));
    }

    BoostPolygon hull;
    boost::geometry::convex_hull(points, hull);
    return Polygon(std::move(hull), false);
}

} // namespace geosets::conversion
