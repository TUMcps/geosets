#pragma once

#include <geosets/sets/polygon.h>
#include <geosets/sets/vpolytope.h>

namespace geosets::conversion {

Polygon vpolytope_to_polygon(const VPolytope& vpolytope);

} // namespace geosets::conversion
