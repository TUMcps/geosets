#pragma once

#include <Eigen/Dense>
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullVertexSet.h>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace geosets::qhull {

// Keeps Qhull plus its backing point buffer alive so libqhull's internal pointers stay valid.
struct ConvexHullHandle {
    std::vector<double> point_storage;
    std::unique_ptr<orgQhull::Qhull> qh_ptr;

    ConvexHullHandle() = default;
    ConvexHullHandle(const ConvexHullHandle&) = delete;
    ConvexHullHandle& operator=(const ConvexHullHandle&) = delete;
    ConvexHullHandle(ConvexHullHandle&&) noexcept = default;
    ConvexHullHandle& operator=(ConvexHullHandle&&) noexcept = default;

    const orgQhull::Qhull& qh() const { return *qh_ptr; }
    orgQhull::Qhull& qh() { return *qh_ptr; }
};

template <typename Derived>
ConvexHullHandle build_convex_hull(const Eigen::MatrixBase<Derived>& points,
                                   const std::string& options = "Qt") {
    using Scalar = typename Derived::Scalar;
    static_assert(std::is_floating_point_v<Scalar>, "Convex hull input must be floating-point");

    const int n = points.rows();
    const int dim = points.cols();

    if (n == 0) throw std::runtime_error("Empty point set");
    if (n < dim + 1) throw std::runtime_error("Not enough points to form a convex hull");

    ConvexHullHandle handle;
    handle.point_storage.reserve(static_cast<size_t>(n) * static_cast<size_t>(dim));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < dim; ++j) {
            handle.point_storage.push_back(static_cast<double>(points(i, j)));
        }
    }

    handle.qh_ptr = std::make_unique<orgQhull::Qhull>();
    handle.qh_ptr->runQhull("", dim, n, handle.point_storage.data(), options.c_str());
    return handle;
}

inline Eigen::MatrixXd vertices_matrix(const ConvexHullHandle& hull) {
    const auto& qh = hull.qh();
    const auto vertices = qh.vertexList();
    const int dim = qh.dimension();
    const int num_vertices = static_cast<int>(vertices.size());

    Eigen::MatrixXd hull_vertices(num_vertices, dim);
    int i = 0;
    for (const auto& vertex : vertices) {
        const auto* p = vertex.point().coordinates();
        for (int d = 0; d < dim; ++d) {
            hull_vertices(i, d) = p[d];
        }
        ++i;
    }

    return hull_vertices;
}

inline double volume(ConvexHullHandle& hull) {
    return hull.qh().volume();
}

} // namespace geosets::qhull
