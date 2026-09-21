#pragma once

#include <geosets/sets/hpolytope.h>
#include <geosets/sets/polygon.h>
#include <geosets/sets/vpolytope.h>

namespace geosets::conversion {

VPolytope hpolytope_to_vpolytope(const HPolytope& hpolytope);
VPolytope polygon_to_vpolytope(const Polygon& polygon);

} // namespace geosets::conversion
