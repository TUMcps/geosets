#pragma once

#include <Eigen/Dense>
#include <set>
#include <string>

namespace geosets::linalg {

std::string vec_to_string(const Eigen::VectorXd& v);

Eigen::MatrixXd vstack(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B);
Eigen::MatrixXd hstack(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B);

size_t rank(const Eigen::MatrixXd& matrix, double tol = 1e-9);

// Helper to convert vertices matrix to a set for order-independent comparison
std::set<std::vector<uint64_t>> vertices_to_set(const Eigen::MatrixXd& vertices);

std::pair<Eigen::MatrixXd, Eigen::VectorXd>
fourier_motzkin_elimination(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, size_t i);

// N-dimensional cross product of (n-1) vectors in R^n.
// Input: (n, n-1) matrix where each column is a vector.
Eigen::VectorXd n_dim_cross_product(const Eigen::MatrixXd& M);

} // namespace geosets::linalg
