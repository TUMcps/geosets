#include <Eigen/Dense>
#include <geosets/core.h>

int main() {
    auto p = geosets::HPolytope::from_unit_box(2);
    p.linear_transform_(Eigen::Matrix2d::Identity() * 1e10);
    std::printf("is_empty: %d\n", p.is_empty());
    auto c = p.center();
    std::printf("center: %f, %f\n", c.x(), c.y());
}
