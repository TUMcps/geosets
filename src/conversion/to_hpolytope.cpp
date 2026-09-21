#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <geosets/conversion/to_hpolytope.h>
#include <geosets/exceptions.h>
#include <geosets/linalg_utils.h>
#include <numeric>
#include <vector>

#include "cdd_utils.h"
#include "geosets/sets/hpolytope.h"

namespace geosets::conversion {

namespace {

void append_inequality(std::vector<Eigen::VectorXd>& rows_A, std::vector<double>& rows_b,
                       const dd_MatrixType* inequalities, dd_rowrange row) {
    const int dim = static_cast<int>(inequalities->colsize - 1);

    Eigen::VectorXd A_row(dim);
    for (int col = 0; col < dim; ++col) {
        A_row(col) = -dd_get_d(inequalities->matrix[row][col + 1]);
    }

    rows_A.push_back(std::move(A_row));
    rows_b.push_back(dd_get_d(inequalities->matrix[row][0]));
}

} // namespace

HPolytope vpolytope_to_hpolytope(const VPolytope& vpolytope) {
    detail::ensure_cdd_initialized();

    const Eigen::MatrixXd& vertices = vpolytope.vertices;
    const int n_vertices = static_cast<int>(vertices.rows());
    const int dim = static_cast<int>(vertices.cols());

    detail::CddMatrixPtr generators = detail::make_cdd_matrix(dd_CreateMatrix(n_vertices, dim + 1));
    if (!generators) {
        throw SetOperationException("Failed to allocate cdd matrix for VPolytope conversion.");
    }
    generators->representation = dd_Generator;

    for (int row = 0; row < n_vertices; ++row) {
        dd_set_si(generators->matrix[row][0], 1);
        for (int col = 0; col < dim; ++col) {
            dd_set_d(generators->matrix[row][col + 1], vertices(row, col));
        }
    }

    dd_ErrorType error = dd_NoError;
    detail::CddPolyhedraPtr polyhedra =
        detail::make_cdd_polyhedra(dd_DDMatrix2Poly(generators.get(), &error));
    if (error != dd_NoError || !polyhedra) {
        detail::throw_cdd_error("Failed VPolytope to HPolytope conversion", error);
    }

    detail::CddMatrixPtr inequalities =
        detail::make_cdd_matrix(dd_CopyInequalities(polyhedra.get()));
    if (!inequalities) {
        throw SetOperationException("Failed to extract inequalities from cdd polyhedron.");
    }

    std::vector<Eigen::VectorXd> rows_A;
    std::vector<double> rows_b;
    rows_A.reserve(static_cast<size_t>(2 * inequalities->rowsize));
    rows_b.reserve(static_cast<size_t>(2 * inequalities->rowsize));

    for (dd_rowrange row = 0; row < inequalities->rowsize; ++row) {
        append_inequality(rows_A, rows_b, inequalities.get(), row);

        if (set_member(row + 1, inequalities->linset)) {
            const std::size_t idx = rows_A.size() - 1;
            rows_A.push_back(-rows_A[idx]);
            rows_b.push_back(-rows_b[idx]);
        }
    }

    if (rows_A.empty()) {
        return HPolytope::make_empty(vpolytope.dim());
    }

    Eigen::MatrixXd A(static_cast<Eigen::Index>(rows_A.size()), dim);
    Eigen::VectorXd b(static_cast<Eigen::Index>(rows_b.size()));
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(rows_A.size()); ++i) {
        A.row(i) = rows_A[static_cast<size_t>(i)].transpose();
        b(i) = rows_b[static_cast<size_t>(i)];
    }

    return HPolytope(std::move(A), std::move(b));
}

HPolytope zonotope_to_hpolytope(const Zonotope& zonotope) {
    const auto n_orig = static_cast<Eigen::Index>(zonotope.dim());
    const auto m = static_cast<Eigen::Index>(zonotope.n_generators());

    // Handle degenerate case via SVD projection
    const bool is_degenerate = zonotope.degenerate();
    Eigen::MatrixXd M_proj;  // n_orig x n_orig orthogonal
    Eigen::VectorXd c_shift; // original center
    Eigen::VectorXd c;
    Eigen::MatrixXd G;
    Eigen::Index n;

    if (is_degenerate) {
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(zonotope.G, Eigen::ComputeFullU);
        const auto r = static_cast<Eigen::Index>(linalg::rank(zonotope.G));

        if (r == 0) {
            // Single point
            Eigen::MatrixXd A(2 * n_orig, n_orig);
            Eigen::VectorXd b(2 * n_orig);
            A.topRows(n_orig) = Eigen::MatrixXd::Identity(n_orig, n_orig);
            A.bottomRows(n_orig) = -Eigen::MatrixXd::Identity(n_orig, n_orig);
            b.head(n_orig) = zonotope.c;
            b.tail(n_orig) = -zonotope.c;
            return HPolytope(std::move(A), std::move(b), false);
        }

        M_proj = svd.matrixU().transpose(); // n_orig x n_orig
        c_shift = zonotope.c;

        // Project: rotate to aligned coordinates, take first r rows
        Eigen::MatrixXd G_rotated = M_proj * zonotope.G; // n_orig x m
        c = Eigen::VectorXd::Zero(r);                    // centered at origin
        G = G_rotated.topRows(r);                        // r x m
        n = r;
    } else {
        c = zonotope.c;
        G = zonotope.G;
        n = n_orig;
    }

    // Number of constraint pairs = C(m, n-1)
    // Compute via combination enumeration
    std::vector<Eigen::VectorXd> normals;
    std::vector<double> offsets_pos;
    std::vector<double> offsets_neg;

    std::vector<Eigen::Index> indices(n - 1);
    std::iota(indices.begin(), indices.end(), Eigen::Index{0});

    while (true) {
        // Build (n, n-1) matrix from selected generator columns
        Eigen::MatrixXd M(n, n - 1);
        for (Eigen::Index j = 0; j < n - 1; ++j) {
            M.col(j) = G.col(indices[j]);
        }

        Eigen::VectorXd cross = linalg::n_dim_cross_product(M);
        double norm = cross.norm();

        if (norm > 1e-12) {
            Eigen::VectorXd normal = cross / norm;
            double delta = (normal.transpose() * G).array().abs().sum();
            double center_dot = normal.dot(c);

            normals.push_back(normal);
            offsets_pos.push_back(center_dot + delta);
            offsets_neg.push_back(-center_dot + delta);
        }

        // Advance to next C(m, n-1) combination
        if (n <= 1) break; // only one (empty) combination
        Eigen::Index i = n - 1;
        while (i > 0) {
            --i;
            ++indices[i];
            if (indices[i] <= m - (n - 1) + i) break;
            if (i == 0) goto done;
        }
        for (Eigen::Index j = i + 1; j < n - 1; ++j) {
            indices[j] = indices[j - 1] + 1;
        }
    }
done:

    const auto h = static_cast<Eigen::Index>(normals.size());
    Eigen::MatrixXd A(2 * h, n);
    Eigen::VectorXd b(2 * h);

    for (Eigen::Index i = 0; i < h; ++i) {
        A.row(i) = normals[i].transpose();
        A.row(i + h) = -normals[i].transpose();
        b(i) = offsets_pos[i];
        b(i + h) = offsets_neg[i];
    }

    // Back-projection for degenerate case
    if (is_degenerate) {
        const Eigen::Index n_extra = n_orig - n;

        Eigen::MatrixXd A_full(2 * h + 2 * n_extra, n_orig);
        Eigen::VectorXd b_full(2 * h + 2 * n_extra);

        // [A, 0; 0, I; 0, -I]
        A_full.topLeftCorner(2 * h, n) = A;
        A_full.topRightCorner(2 * h, n_extra).setZero();
        A_full.block(2 * h, 0, n_extra, n).setZero();
        A_full.block(2 * h, n, n_extra, n_extra).setIdentity();
        A_full.block(2 * h + n_extra, 0, n_extra, n).setZero();
        A_full.block(2 * h + n_extra, n, n_extra, n_extra) =
            -Eigen::MatrixXd::Identity(n_extra, n_extra);

        b_full.head(2 * h) = b;
        b_full.tail(2 * n_extra).setZero();

        // Rotate back: A_full in rotated coords → original coords via M_proj^T = U
        A_full = A_full * M_proj;

        // Re-center
        b_full += A_full * c_shift;

        return HPolytope(std::move(A_full), std::move(b_full), false);
    }

    return HPolytope(std::move(A), std::move(b), false);
}

} // namespace geosets::conversion
