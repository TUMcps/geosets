#pragma once

#include <Eigen/Dense>

namespace geosets {

struct Support {
    double value;
    Eigen::VectorXd vector;
};

struct Hyperplane {
    // Represents the hyperplane {x | a^T x = b}
    Eigen::VectorXd a;
    double b;
};

} // namespace geosets
