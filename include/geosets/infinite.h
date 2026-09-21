#pragma once

#include "geosets/exceptions.h"
#include <optional>

namespace geosets {

/**
 * @brief Checks if one of the sets is infinite and returns the other set as the result if so.
 * Required for e.g. intersection.
 */
template <typename SetType>
inline std::optional<SetType> check_infinite_identity(const SetType& set1, const SetType& set2) {
    if (set1.is_infinite()) return set2.copy();
    if (set2.is_infinite()) return set1.copy();
    return std::nullopt;
}

/**
 * @brief Checks if one of the sets is infinite and returns an infinite set as the result if so.
 * Required for e.g. minkowski sum.
 */
template <typename SetType>
inline std::optional<SetType> check_infinite_dominant(const SetType& set1, const SetType& set2) {
    if (set1.is_infinite() || set2.is_infinite()) {
        return SetType::make_infinite(set1.dim());
    }
    return std::nullopt;
}

/**
 * @brief Checks if one of the sets is infinite for subtractive operations.
 * Required for e.g. set difference.
 */
template <typename SetType>
inline std::optional<SetType> check_infinite_subtractive(const SetType& set1, const SetType& set2) {
    if (set2.is_infinite()) return SetType::make_empty(set1.dim());
    if (set1.is_infinite()) {
        throw SetOperationException(
            "set_difference: Subtracting from an infinite set is not supported");
    }
    return std::nullopt;
}

} // namespace geosets
