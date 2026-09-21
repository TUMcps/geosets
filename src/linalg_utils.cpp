#include <cassert>
#include <cstddef>
#include <geosets/linalg_utils.h>

namespace geosets::linalg {

std::string vec_to_string(const Eigen::VectorXd& v) {
    std::stringstream ss;
    Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n", "[", "]");
    ss << "Vec:\n" << v.format(fmt);
    return ss.str();
}

Eigen::MatrixXd vstack(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
    assert(A.cols() == B.cols());
    Eigen::MatrixXd result(A.rows() + B.rows(), A.cols());
    result << A, B;
    return result;
}

Eigen::MatrixXd hstack(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
    assert(A.rows() == B.rows());
    Eigen::MatrixXd result(A.rows(), A.cols() + B.cols());
    result << A, B;
    return result;
}

size_t rank(const Eigen::MatrixXd& matrix, double tol) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix);
    const auto& singular_values = svd.singularValues();

    int rank = 0;
    for (int i = 0; i < singular_values.size(); ++i) {
        if (singular_values(i) > tol) {
            rank++;
        }
    }
    return rank;
}

std::pair<Eigen::MatrixXd, Eigen::VectorXd>
fourier_motzkin_elimination(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, size_t i) {
    // For stability, set all values closer than a given tolerance to exactly 0
    const double atol = 1e-8;
    Eigen::MatrixXd A_proj = A;
    A_proj = (A_proj.array().abs() < atol).select(0.0, A_proj);

    // Divide the i-th column into positive, zero, and negative entries
    std::vector<size_t> Z, P, N;
    for (size_t row = 0; row < static_cast<size_t>(A_proj.rows()); ++row) {
        if (A_proj(row, i) == 0) {
            Z.push_back(row);
        } else if (A_proj(row, i) > 0) {
            P.push_back(row);
        } else {
            N.push_back(row);
        }
    }

    const size_t m = b.size();
    const size_t n_rows = Z.size() + P.size() * N.size();
    Eigen::MatrixXd U = Eigen::MatrixXd::Zero(n_rows, m);

    // Deal with Z - these constraints remain unchanged
    for (size_t j = 0; j < Z.size(); ++j) {
        U(j, Z[j]) = 1.0;
    }

    // Deal with N x P - Cartesian product
    size_t idx = Z.size();
    for (size_t n_idx : N) {
        for (size_t p_idx : P) {
            U(idx, n_idx) = A_proj(p_idx, i);
            U(idx, p_idx) = -A_proj(n_idx, i);
            ++idx;
        }
    }
    A_proj = U * A_proj;

    // Remove column i
    Eigen::MatrixXd A_result(A_proj.rows(), A_proj.cols() - 1);
    if (i > 0) {
        A_result.leftCols(i) = A_proj.leftCols(i);
    }
    if (i < static_cast<size_t>(A_proj.cols()) - 1) {
        A_result.rightCols(A_proj.cols() - i - 1) = A_proj.rightCols(A_proj.cols() - i - 1);
    }

    Eigen::VectorXd b_proj = U * b;

    return {A_result, b_proj};
}
Eigen::VectorXd n_dim_cross_product(const Eigen::MatrixXd& M) {
    assert(M.cols() == M.rows() - 1);
    const Eigen::Index n = M.rows();
    Eigen::VectorXd result(n);

    for (Eigen::Index i = 0; i < n; ++i) {
        // Build (n-1)x(n-1) submatrix by removing row i
        Eigen::MatrixXd sub(n - 1, n - 1);
        Eigen::Index row = 0;
        for (Eigen::Index j = 0; j < n; ++j) {
            if (j == i) continue;
            sub.row(row++) = M.row(j);
        }
        // (-1)^(i+2) = (-1)^i since (-1)^2 = 1
        result(i) = ((i % 2 == 0) ? 1.0 : -1.0) * sub.determinant();
    }

    return result;
}

} // namespace geosets::linalg
