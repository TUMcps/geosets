#pragma once

#include <Eigen/Dense>
#include <cstddef>
#include <geosets/parallel.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/map.h>

// Sets are passed by pointer so the original Python-owned C++ objects are
// reused across calls. This preserves any cached state (e.g. Chebyshev center
// on HPolytope) instead of operating on per-call copies.

#define GS_BIND_BATCHED_VOID(NAME)                                                                 \
    m.def(("batched_" #NAME "_" + typestr).c_str(),                                                \
          [](const std::vector<SetType*>& sets) {                                                  \
              return geosets::parallel_execute(sets,                                               \
                                               [](SetType* s, size_t) { return s->NAME(); });      \
          },                                                                                       \
          nb::arg("sets"))

template <typename SetType>
void parallelized_bindings(nanobind::module_& m, const std::string& typestr) {
    namespace nb = nanobind;

    GS_BIND_BATCHED_VOID(volume);
    GS_BIND_BATCHED_VOID(is_empty);
    GS_BIND_BATCHED_VOID(bounded);
    GS_BIND_BATCHED_VOID(degenerate);
    GS_BIND_BATCHED_VOID(to_vertices);

    m.def(("batched_center_" + typestr).c_str(),
          [](const std::vector<SetType*>& sets, bool validate) {
              return geosets::parallel_execute(
                  sets, [&](SetType* s, size_t) { return s->center(validate); });
          },
          nb::arg("sets"), nb::arg("validate") = false);

    // Per-set direction: directions.row(i).transpose() is the query for sets[i].
    m.def(("batched_support_function_" + typestr).c_str(),
          [](const std::vector<SetType*>& sets, const Eigen::MatrixXd& directions, bool validate) {
              return geosets::parallel_execute(sets, [&](SetType* s, size_t i) {
                  return s->support_function(directions.row(i).transpose(), validate);
              });
          },
          nb::arg("sets"), nb::arg("directions"), nb::arg("validate") = false);

    m.def(("batched_tangent_hyperplane_" + typestr).c_str(),
          [](const std::vector<SetType*>& sets, const Eigen::MatrixXd& directions, bool validate) {
              return geosets::parallel_execute(sets, [&](SetType* s, size_t i) {
                  return s->tangent_hyperplane(directions.row(i).transpose(), validate);
              });
          },
          nb::arg("sets"), nb::arg("directions"), nb::arg("validate") = false);

    m.def(("batched_boundary_point_" + typestr).c_str(),
          [](const std::vector<SetType*>& sets, const Eigen::MatrixXd& directions, bool validate) {
              return geosets::parallel_execute(sets, [&](SetType* s, size_t i) {
                  return s->boundary_point(directions.row(i).transpose(), validate);
              });
          },
          nb::arg("sets"), nb::arg("directions"), nb::arg("validate") = false);

    // Per-set point: points.row(i) is the query for sets[i] (shape: batch x dim).
    m.def(("batched_contains_point_" + typestr).c_str(),
          [](const std::vector<SetType*>& sets, const Eigen::MatrixXd& points, double tol,
             bool validate) {
              return geosets::parallel_execute(sets, [&](SetType* s, size_t i) {
                  return s->contains_point(points.row(i).transpose(), tol, validate);
              });
          },
          nb::arg("sets"), nb::arg("points"), nb::arg("tol") = 1e-9, nb::arg("validate") = false);
}

#undef GS_BIND_BATCHED_VOID
