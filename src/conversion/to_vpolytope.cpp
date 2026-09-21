#include <Eigen/Dense>
#include <cmath>
#include <geosets/conversion/to_vpolytope.h>
#include <geosets/exceptions.h>
#include <vector>

#include "cdd_utils.h"

namespace geosets::conversion {

namespace {

constexpr double kTol = 1e-9;

} // namespace

VPolytope hpolytope_to_vpolytope(const HPolytope& hpolytope) {
    detail::ensure_cdd_initialized();

    const Eigen::MatrixXd& A = hpolytope.A;
    const Eigen::VectorXd& b = hpolytope.b;

    const int m = static_cast<int>(A.rows());
    const int n = static_cast<int>(A.cols());

    // Removed checks for performance
    // if (hpolytope.is_empty()) {
    //     throw SetOperationException("Cannot convert empty HPolytope to VPolytope.");
    // }

    // if (!hpolytope.bounded()) {
    //     throw SetOperationException(
    //         "HPolytope to VPolytope conversion requires a bounded polytope.");
    // }

    detail::CddMatrixPtr inequalities = detail::make_cdd_matrix(dd_CreateMatrix(m, n + 1));
    if (!inequalities) {
        throw SetOperationException("Failed to allocate cdd matrix for HPolytope conversion.");
    }
    inequalities->representation = dd_Inequality;

    for (int row = 0; row < m; ++row) {
        dd_set_d(inequalities->matrix[row][0], b(row));
        for (int col = 0; col < n; ++col) {
            dd_set_d(inequalities->matrix[row][col + 1], -A(row, col));
        }
    }

    dd_ErrorType error = dd_NoError;
    detail::CddPolyhedraPtr polyhedra =
        detail::make_cdd_polyhedra(dd_DDMatrix2Poly(inequalities.get(), &error));
    if (error != dd_NoError || !polyhedra) {
        detail::throw_cdd_error("Failed HPolytope to VPolytope conversion", error);
    }

    detail::CddMatrixPtr generators = detail::make_cdd_matrix(dd_CopyGenerators(polyhedra.get()));
    if (!generators) {
        throw SetOperationException("Failed to extract generators from cdd polyhedron.");
    }

    std::vector<Eigen::VectorXd> vertices;
    vertices.reserve(static_cast<size_t>(generators->rowsize));

    for (dd_rowrange row = 0; row < generators->rowsize; ++row) {
        if (set_member(row + 1, generators->linset)) {
            throw SetOperationException("Bounded HPolytope conversion produced line generators.");
        }

        const double homogeneous = dd_get_d(generators->matrix[row][0]);
        if (std::abs(homogeneous) <= kTol) {
            throw SetOperationException("Bounded HPolytope conversion produced ray generators.");
        }

        Eigen::VectorXd vertex(n);
        for (int col = 0; col < n; ++col) {
            const double value = dd_get_d(generators->matrix[row][col + 1]);
            vertex(col) = value / homogeneous;
        }
        vertices.push_back(std::move(vertex));
    }

    if (vertices.empty()) {
        return VPolytope::make_empty(static_cast<std::size_t>(n));
    }

    Eigen::MatrixXd vertices_matrix(static_cast<Eigen::Index>(vertices.size()), n);
    for (Eigen::Index i = 0; i < vertices_matrix.rows(); ++i) {
        vertices_matrix.row(i) = vertices[static_cast<size_t>(i)].transpose();
    }

    return VPolytope(std::move(vertices_matrix), false);
}

VPolytope polygon_to_vpolytope(const Polygon& polygon) {
    return VPolytope(polygon.to_vertices(), false);
}

} // namespace geosets::conversion
