#pragma once

// The files need to be included exactly in this order!
// clang-format off
#include <cddlib/setoper.h>
#include <cddlib/cdd.h>
// clang-format on
#include <cstdlib>
#include <geosets/exceptions.h>
#include <memory>
#include <mutex>
#include <string>

namespace geosets::conversion::detail {

inline void ensure_cdd_initialized() {
    static std::once_flag init_flag;
    std::call_once(init_flag, []() {
        dd_set_global_constants();
        std::atexit(dd_free_global_constants);
    });
}

struct CddMatrixDeleter {
    void operator()(dd_MatrixType* matrix) const {
        if (matrix != nullptr) {
            dd_FreeMatrix(matrix);
        }
    }
};

struct CddPolyhedraDeleter {
    void operator()(dd_PolyhedraType* poly) const {
        if (poly != nullptr) {
            dd_FreePolyhedra(poly);
        }
    }
};

using CddMatrixPtr = std::unique_ptr<dd_MatrixType, CddMatrixDeleter>;
using CddPolyhedraPtr = std::unique_ptr<dd_PolyhedraType, CddPolyhedraDeleter>;

inline CddMatrixPtr make_cdd_matrix(dd_MatrixPtr matrix) {
    return CddMatrixPtr(matrix);
}

inline CddPolyhedraPtr make_cdd_polyhedra(dd_PolyhedraPtr poly) {
    return CddPolyhedraPtr(poly);
}

[[noreturn]] inline void throw_cdd_error(const std::string& context, dd_ErrorType error) {
    throw SetOperationException(
        context + " (cdd error code: " + std::to_string(static_cast<int>(error)) + ").");
}

} // namespace geosets::conversion::detail
