#pragma once

#include <geosets/conversion/to_hpolytope.h>
#include <geosets/conversion/to_polygon.h>
#include <geosets/conversion/to_vpolytope.h>
#include <geosets/parallel.h>
#include <geosets/sampling/gaussian_combined.h>
#include <geosets/sampling/gaussian_rdhr.h>
#include <geosets/sampling/gaussian_rejection.h>
#include <geosets/sets/hpolytope.h>

#include <ctime>
#include <geosets/sets/interface.h>
#include <geosets/sets/interval.h>
#include <geosets/sets/polygon.h>
#include <geosets/sets/vpolytope.h>
#include <geosets/sets/zonotope.h>

namespace geosets {

inline void random_seed(unsigned int seed = 0) {
    if (seed == 0) seed = static_cast<unsigned int>(time(0));
    seed_rng(seed);
}

inline void check_sets_dim(const GeometricSetInterface* s1, const GeometricSetInterface* s2) {
    // Note: This could potentially be slightly slower than checking directly due to function call.
    // But very likely negligible due to compiler otpimizations
    if (s1->dim() != s2->dim()) {
        throw DimensionMismatchException(s1->dim(), s2->dim());
    }
}

} // namespace geosets
