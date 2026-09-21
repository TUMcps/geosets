#include <Eigen/Dense>
#include <set>

// Helper to convert vertices matrix to a set for order-independent comparison

inline std::set<std::vector<uint64_t>> vertices_to_set(const Eigen::MatrixXd& vertices) {
    std::set<std::vector<uint64_t>> result;
    for (int i = 0; i < vertices.rows(); ++i) {
        std::vector<uint64_t> row;
        for (int j = 0; j < vertices.cols(); ++j) {
            uint64_t bits;
            double val = vertices(i, j);
            std::memcpy(&bits, &val, sizeof(bits));
            row.push_back(bits);
        }
        result.insert(row);
    }
    return result;
}
