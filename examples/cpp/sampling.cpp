#include "geosets/sampling/gaussian_rdhr.h"
#include <cstddef>
#include <geosets/core.h>

int main() {
    auto offset = Eigen::VectorXd::Ones(2) * 6;
    auto p = geosets::HPolytope::from_unit_box(2).translate(offset);
    auto loc = Eigen::VectorXd::Zero(2);
    auto scale = Eigen::VectorXd::Ones(2);

    // Note: parameter name should be n_samples instead of samples_per_set
    auto samples = geosets::sampling::gaussian_rdhr_samples(p, loc, scale, 10);

    // print samples
    for (size_t i = 0; i < samples.rows(); ++i) {
        for (size_t j = 0; j < samples.cols(); ++j) {
            std::cout << samples(i, j) << " ";
        }
        std::cout << std::endl;
    }
}
