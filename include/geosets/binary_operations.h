#pragma once

#include "geosets/core.h"
#include <geosets/infinite.h>

namespace geosets {

// TODO: The check_inifinite_ functions could be just moved to the function body, since they are
// only used here anyway.
//  This might make it easier to read.

// Generic intersection predicate for any set type
template <typename SetType> bool do_intersect(const SetType& s1, const SetType& s2) {
    geosets::check_sets_dim(&s1, &s2);
    if (s1.is_empty() || s2.is_empty()) return false;
    if (s1.is_infinite() || s2.is_infinite()) return true;
    return do_intersect_impl(s1, s2);
}

// Generic intersection for any set type
template <typename SetType> SetType intersection(const SetType& s1, const SetType& s2) {
    geosets::check_sets_dim(&s1, &s2);
    if (auto result = check_infinite_identity(s1, s2)) {
        return result.value();
    }
    return intersection_impl(s1, s2);
}

// Generic minkowski_sum for any set type
template <typename SetType> SetType minkowski_sum(const SetType& s1, const SetType& s2) {
    geosets::check_sets_dim(&s1, &s2);
    if (auto result = check_infinite_dominant(s1, s2)) {
        return result.value();
    }
    return minkowski_sum_impl(s1, s2);
}

// Generic set difference for any set type
template <typename SetType> SetType set_difference(const SetType& s1, const SetType& s2) {
    geosets::check_sets_dim(&s1, &s2);
    if (auto result = check_infinite_subtractive(s1, s2)) {
        return result.value();
    }
    if (s2.is_empty()) return s1.copy();
    return set_difference_impl(s1, s2);
}

} // namespace geosets
