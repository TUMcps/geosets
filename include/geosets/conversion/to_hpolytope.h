#pragma once

#include <geosets/sets/hpolytope.h>
#include <geosets/sets/vpolytope.h>
#include <geosets/sets/zonotope.h>

namespace geosets::conversion {

HPolytope vpolytope_to_hpolytope(const VPolytope& vpolytope);
HPolytope zonotope_to_hpolytope(const Zonotope& zonotope);

} // namespace geosets::conversion
