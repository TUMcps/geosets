#include "geosets/binary_operations.h"
#include "geosets/core.h"
#include "geosets/exceptions.h"
#include "geosets/structs.h"
#include "parallelized_bindings.h"
#include "sampling_bindings.h"
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <cstddef>
#include <geosets/sets/hpolytope.h>
#include <geosets/sets/interface.h>
#include <geosets/sets/interval.h>
#include <geosets/sets/polygon.h>
#include <geosets/sets/vpolytope.h>
#include <geosets/sets/zonotope.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/unique_ptr.h>
#include <nanobind/stl/vector.h>
#include <tuple>

namespace nb = nanobind;

// _py_core must be the name of the nanobind module defined in the
// CMakeLists.txt file
NB_MODULE(_gs_py_core, m) {
    m.doc() = "Python bindings for the geosets library.";

    // -- General bindings --
    nb::class_<geosets::Support>(m, "Support",
                                 "Struct representing the result of a support function evaluation.")
        .def_ro("value", &geosets::Support::value)
        .def_ro("vector", &geosets::Support::vector, nb::rv_policy::reference_internal);

    nb::class_<geosets::Hyperplane>(m, "Hyperplane", "Struct representing hyperplane a^T x = b.")
        .def_ro("b", &geosets::Hyperplane::b)
        .def_ro("a", &geosets::Hyperplane::a, nb::rv_policy::reference_internal);

    // -- Exception bindings --
    nb::exception<geosets::NotImplemented> _1(m, "NotImplemented");
    nb::exception<geosets::SetOperationException> _2(m, "SetOperationException");
    nb::exception<geosets::DimensionMismatchException> _3(m, "DimensionMismatchException");
    nb::exception<geosets::FailedOptimizationProblem> _4(m, "FailedOptimizationProblem");
    nb::exception<geosets::InfeasibleOptimizationProblem> _5(m, "InfeasibleOptimizationProblem");
    nb::exception<geosets::NotOptimalOptimizationProblem> _6(m, "NotOptimalOptimizationProblem");
    nb::exception<geosets::DegenerateSetException> _7(m, "DegenerateSetException");
    nb::exception<geosets::InfiniteSetException> _8(m, "InfiniteSetException");
    nb::exception<geosets::EmptySetException> _9(m, "EmptySetException");

    // -- Set implementation bindings --
    nb::class_<geosets::GeometricSetInterface>(m, "GeometricSet", "Geometric set interface class.")
        .def("__str__", &geosets::GeometricSetInterface::str)
        .def("dim", &geosets::GeometricSetInterface::dim)
        .def("is_empty", &geosets::GeometricSetInterface::is_empty)
        .def("bounded", &geosets::GeometricSetInterface::bounded)
        .def("degenerate", &geosets::GeometricSetInterface::degenerate)
        .def("to_vertices", &geosets::GeometricSetInterface::to_vertices)
        .def("center", &geosets::GeometricSetInterface::center, nb::arg("validate") = true)
        .def("contains_point", &geosets::GeometricSetInterface::contains_point, nb::arg("point"),
             nb::arg("tol") = 1e-9, nb::arg("validate") = true)
        .def("volume", &geosets::GeometricSetInterface::volume)
        .def("support_function", &geosets::GeometricSetInterface::support_function,
             nb::arg("direction"), nb::arg("validate") = true)
        .def("tangent_hyperplane", &geosets::GeometricSetInterface::tangent_hyperplane,
             nb::arg("direction"), nb::arg("validate") = true)
        .def("boundary_point", &geosets::GeometricSetInterface::boundary_point,
             nb::arg("direction"), nb::arg("validate") = true)
        .def("line_intersection", &geosets::GeometricSetInterface::line_intersection,
             nb::arg("point"), nb::arg("direction"), nb::arg("validate") = true)
        .def("translate_", &geosets::GeometricSetInterface::translate_, nb::arg("vector"))
        .def("linear_transform_", &geosets::GeometricSetInterface::linear_transform_,
             nb::arg("matrix"))
        .def(
            "plot",
            [](geosets::GeometricSetInterface& self, nb::object ax, std::tuple<int, int> axis,
               bool show, double fill_opacity, nb::kwargs kwargs) {
                throw std::runtime_error(
                    "plot() is implemented in Python. Import _py_core to attach it.");
            },
            nb::arg("ax") = nb::none(), nb::arg("axis") = nb::make_tuple(0, 1),
            nb::arg("show") = false, nb::arg("fill_opacity") = 0.0,
            nb::arg("kwargs") = nb::kwargs())
        // -- These definitions are just here to generate the python stubs for the GeometricSet base
        // class -- Maybe there is a better method to do this?
        .def("copy",
             [](geosets::GeometricSetInterface& self) -> geosets::GeometricSetInterface* {
                 throw nb::type_error("Method 'copy' is not implemented in the base interface.");
             })
        .def(
            "translate",
            [](geosets::GeometricSetInterface& self,
               const Eigen::VectorXd& vector) -> geosets::GeometricSetInterface* {
                throw nb::type_error(
                    "Method 'translate' is not implemented in the base interface.");
            },
            nb::arg("vector"))
        .def(
            "linear_transform",
            [](geosets::GeometricSetInterface& self,
               const Eigen::MatrixXd& matrix) -> geosets::GeometricSetInterface* {
                throw nb::type_error(
                    "Method 'linear_transform' is not implemented in the base interface.");
            },
            nb::arg("vector"));

    // Note: methods that return an instance of the class, like translate, or linear_transform need
    // to be bound for the sub-classes individually!

    nb::class_<geosets::Interval, geosets::GeometricSetInterface>(m, "Interval")
        .def(nb::init<Eigen::VectorXd, Eigen::VectorXd, bool>(), nb::arg("lb"), nb::arg("ub"),
             nb::arg("validate") = true)
        .def("__getstate__",
             [](const geosets::Interval& self) { return std::make_tuple(self.lb, self.ub); })
        .def("__setstate__",
             [](geosets::Interval& self, std::tuple<Eigen::VectorXd, Eigen::VectorXd> state) {
                 const auto& [lb, ub] = state;
                 new (&self) geosets::Interval(lb, ub);
             })
        .def_ro("lb", &geosets::Interval::lb, nb::rv_policy::reference_internal)
        .def_ro("ub", &geosets::Interval::ub, nb::rv_policy::reference_internal)
        .def_static("from_unit_box", &geosets::Interval::from_unit_box, nb::arg("dim"))
        .def_static("from_random", &geosets::Interval::from_random, nb::arg("dim"))
        .def("copy", &geosets::Interval::copy_to_ptr, nb::rv_policy::take_ownership)
        .def("translate", &geosets::Interval::translate_ptr, nb::arg("vector"),
             nb::rv_policy::take_ownership)
        .def("linear_transform", &geosets::Interval::linear_transform_ptr, nb::arg("matrix"),
             nb::rv_policy::take_ownership);

    // Binary operations for Interval
    // They need to be static_cast due to the overloaded implementations
    m.def("intersection",
          static_cast<geosets::Interval (*)(const geosets::Interval&, const geosets::Interval&)>(
              &geosets::intersection),
          nb::arg("s1"), nb::arg("s2"));
    m.def("do_intersect",
          static_cast<bool (*)(const geosets::Interval&, const geosets::Interval&)>(
              &geosets::do_intersect),
          nb::arg("s1"), nb::arg("s2"));
    m.def("minkowski_sum",
          static_cast<geosets::Interval (*)(const geosets::Interval&, const geosets::Interval&)>(
              &geosets::minkowski_sum),
          nb::arg("s1"), nb::arg("s2"));

    nb::class_<geosets::HPolytope, geosets::GeometricSetInterface>(m, "HPolytope")
        .def(nb::init<Eigen::MatrixXd, Eigen::VectorXd, bool>(), nb::arg("A"), nb::arg("b"),
             nb::arg("validate") = true)
        .def("__getstate__",
             [](const geosets::HPolytope& self) { return std::make_tuple(self.A, self.b); })
        .def("__setstate__",
             [](geosets::HPolytope& self, std::tuple<Eigen::MatrixXd, Eigen::VectorXd> state) {
                 const auto& [A, b] = state;
                 new (&self) geosets::HPolytope(A, b);
             })
        // We need custom copy and deepcopy to copy the cached center
        .def("__copy__", [](const geosets::HPolytope& self) { return geosets::HPolytope(self); })
        .def(
            "__deepcopy__",
            [](const geosets::HPolytope& self, nb::handle /*memo*/) {
                return geosets::HPolytope(self);
            },
            nb::arg("memo"))
        .def_ro("A", &geosets::HPolytope::A, nb::rv_policy::reference_internal)
        .def_ro("b", &geosets::HPolytope::b, nb::rv_policy::reference_internal)
        .def_static("from_unit_box", &geosets::HPolytope::from_unit_box, nb::arg("dim"))
        .def_static("from_random", &geosets::HPolytope::from_random, nb::arg("dim"),
                    nb::arg("n_constraints"))
        .def("copy", &geosets::HPolytope::copy_to_ptr, nb::rv_policy::take_ownership)
        .def("translate", &geosets::HPolytope::translate_ptr, nb::arg("vector"),
             nb::rv_policy::take_ownership)
        .def("linear_transform", &geosets::HPolytope::linear_transform_ptr, nb::arg("matrix"),
             nb::rv_policy::take_ownership);

    // Binary operations for HPolytope
    m.def("intersection",
          static_cast<geosets::HPolytope (*)(const geosets::HPolytope&, const geosets::HPolytope&)>(
              &geosets::intersection),
          nb::arg("s1"), nb::arg("s2"));
    m.def("do_intersect",
          static_cast<bool (*)(const geosets::HPolytope&, const geosets::HPolytope&)>(
              &geosets::do_intersect),
          nb::arg("s1"), nb::arg("s2"));
    m.def("minkowski_sum",
          static_cast<geosets::HPolytope (*)(const geosets::HPolytope&, const geosets::HPolytope&)>(
              &geosets::minkowski_sum),
          nb::arg("s1"), nb::arg("s2"));

    nb::class_<geosets::VPolytope, geosets::GeometricSetInterface>(m, "VPolytope")
        .def(nb::init<Eigen::MatrixXd, bool>(), nb::arg("vertices"), nb::arg("validate") = true)
        .def("__getstate__", [](const geosets::VPolytope& self) { return self.to_vertices(); })
        .def("__setstate__",
             [](geosets::VPolytope& self, Eigen::MatrixXd vertices) {
                 new (&self) geosets::VPolytope(vertices);
             })
        .def_ro("vertices", &geosets::VPolytope::vertices, nb::rv_policy::reference_internal)
        .def_static("from_unit_box", &geosets::VPolytope::from_unit_box, nb::arg("dim"))
        .def_static("from_random", &geosets::VPolytope::from_random, nb::arg("dim"),
                    nb::arg("n_vertices"))
        .def("copy", &geosets::VPolytope::copy_to_ptr, nb::rv_policy::take_ownership)
        .def("translate", &geosets::VPolytope::translate_ptr, nb::arg("vector"),
             nb::rv_policy::take_ownership)
        .def("linear_transform", &geosets::VPolytope::linear_transform_ptr, nb::arg("matrix"),
             nb::rv_policy::take_ownership);

    // Binary operations for VPolytope
    m.def("intersection",
          static_cast<geosets::VPolytope (*)(const geosets::VPolytope&, const geosets::VPolytope&)>(
              &geosets::intersection),
          nb::arg("s1"), nb::arg("s2"));
    m.def("do_intersect",
          static_cast<bool (*)(const geosets::VPolytope&, const geosets::VPolytope&)>(
              &geosets::do_intersect),
          nb::arg("s1"), nb::arg("s2"));
    m.def("minkowski_sum",
          static_cast<geosets::VPolytope (*)(const geosets::VPolytope&, const geosets::VPolytope&)>(
              &geosets::minkowski_sum),
          nb::arg("s1"), nb::arg("s2"));

    nb::class_<geosets::Zonotope, geosets::GeometricSetInterface>(m, "Zonotope")
        .def(nb::init<Eigen::VectorXd, Eigen::MatrixXd, bool>(), nb::arg("c"), nb::arg("G"),
             nb::arg("validate") = true)
        .def("__getstate__",
             [](const geosets::Zonotope& self) { return std::make_tuple(self.c, self.G); })
        .def("__setstate__",
             [](geosets::Zonotope& self, std::tuple<Eigen::VectorXd, Eigen::MatrixXd> state) {
                 const auto& [c, G] = state;
                 new (&self) geosets::Zonotope(c, G);
             })
        .def_ro("c", &geosets::Zonotope::c, nb::rv_policy::reference_internal)
        .def_ro("G", &geosets::Zonotope::G, nb::rv_policy::reference_internal)
        .def_static("from_unit_box", &geosets::Zonotope::from_unit_box, nb::arg("dim"))
        .def_static("from_random", &geosets::Zonotope::from_random, nb::arg("dim"),
                    nb::arg("n_generators"))
        .def("n_generators", &geosets::Zonotope::n_generators)
        .def("order", &geosets::Zonotope::order)
        .def("zonotope_norm", &geosets::Zonotope::zonotope_norm, nb::arg("point"))
        .def("copy", &geosets::Zonotope::copy_to_ptr, nb::rv_policy::take_ownership)
        .def("translate", &geosets::Zonotope::translate_ptr, nb::arg("vector"),
             nb::rv_policy::take_ownership)
        .def("linear_transform", &geosets::Zonotope::linear_transform_ptr, nb::arg("matrix"),
             nb::rv_policy::take_ownership);

    // Binary operations for Zonotope
    m.def("intersection",
          static_cast<geosets::Zonotope (*)(const geosets::Zonotope&, const geosets::Zonotope&)>(
              &geosets::intersection),
          nb::arg("s1"), nb::arg("s2"));
    m.def("do_intersect",
          static_cast<bool (*)(const geosets::Zonotope&, const geosets::Zonotope&)>(
              &geosets::do_intersect),
          nb::arg("s1"), nb::arg("s2"));
    m.def("minkowski_sum",
          static_cast<geosets::Zonotope (*)(const geosets::Zonotope&, const geosets::Zonotope&)>(
              &geosets::minkowski_sum),
          nb::arg("s1"), nb::arg("s2"));

    nb::class_<geosets::BoostPolygon> boost_polygon_view(m, "_BoostPolygon");
    nb::class_<geosets::Polygon, geosets::GeometricSetInterface>(m, "Polygon")
        .def(nb::init<Eigen::MatrixXd, bool>(), nb::arg("vertices"), nb::arg("validate") = true)
        .def("__getstate__", [](const geosets::Polygon& self) { return self.to_vertices(); })
        .def("__setstate__",
             [](geosets::Polygon& self, Eigen::MatrixXd vertices) {
                 new (&self) geosets::Polygon(vertices);
             })
        .def_ro("boost_polygon", &geosets::Polygon::boost_polygon,
                nb::rv_policy::reference_internal)
        .def_static("from_unit_box", &geosets::Polygon::from_unit_box, nb::arg("dim"))
        .def_static("from_random", &geosets::Polygon::from_random, nb::arg("dim"),
                    nb::arg("n_vertices"))
        .def("copy", &geosets::Polygon::copy_to_ptr, nb::rv_policy::take_ownership)
        .def("translate", &geosets::Polygon::translate_ptr, nb::arg("vector"),
             nb::rv_policy::take_ownership)
        .def("linear_transform", &geosets::Polygon::linear_transform_ptr, nb::arg("matrix"),
             nb::rv_policy::take_ownership);

    // Binary operations for Polygon
    m.def("intersection",
          static_cast<geosets::Polygon (*)(const geosets::Polygon&, const geosets::Polygon&)>(
              &geosets::intersection),
          nb::arg("s1"), nb::arg("s2"));
    m.def("do_intersect",
          static_cast<bool (*)(const geosets::Polygon&, const geosets::Polygon&)>(
              &geosets::do_intersect),
          nb::arg("s1"), nb::arg("s2"));
    m.def("minkowski_sum",
          static_cast<geosets::Polygon (*)(const geosets::Polygon&, const geosets::Polygon&)>(
              &geosets::minkowski_sum),
          nb::arg("s1"), nb::arg("s2"));
    m.def("set_difference",
          static_cast<geosets::Polygon (*)(const geosets::Polygon&, const geosets::Polygon&)>(
              &geosets::set_difference),
          nb::arg("s1"), nb::arg("s2"));

    // --- Conversion bindings ---
    m.def("hpolytope_to_vpolytope", &geosets::conversion::hpolytope_to_vpolytope,
          nb::arg("hpolytope"));
    m.def("vpolytope_to_hpolytope", &geosets::conversion::vpolytope_to_hpolytope,
          nb::arg("vpolytope"));
    m.def("zonotope_to_hpolytope", &geosets::conversion::zonotope_to_hpolytope,
          nb::arg("zonotope"));
    m.def("polygon_to_vpolytope", &geosets::conversion::polygon_to_vpolytope, nb::arg("polygon"));
    m.def("vpolytope_to_polygon", &geosets::conversion::vpolytope_to_polygon, nb::arg("vpolytope"));

    // --- Function bindings ---
    m.def("random_seed", &geosets::random_seed, nb::arg("seed") = 0);

    // --- Batched computation bindings ---
    // The template functions need to be bound for each set (Python does not have templates)
    m.def("set_omp_num_threads", &geosets::set_omp_num_threads, nb::arg("nthreads"));
    parallelized_bindings<geosets::HPolytope>(m, "hpolytope");
    parallelized_bindings<geosets::VPolytope>(m, "vpolytope");
    parallelized_bindings<geosets::Polygon>(m, "polygon");
    parallelized_bindings<geosets::Zonotope>(m, "zonotope");
    parallelized_bindings<geosets::Interval>(m, "interval");

    nb::module_ sampling = m.def_submodule("sampling", "Gaussian sampling routines.");
    sampling_bindings<geosets::HPolytope>(sampling, "hpolytope");
    sampling_bindings<geosets::VPolytope>(sampling, "vpolytope");
    sampling_bindings<geosets::Polygon>(sampling, "polygon");
    sampling_bindings<geosets::Zonotope>(sampling, "zonotope");
    sampling_bindings<geosets::Interval>(sampling, "interval");
}
