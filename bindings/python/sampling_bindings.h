#pragma once

#include <geosets/exceptions.h>
#include <geosets/parallel.h>
#include <geosets/sampling/common.h>
#include <geosets/sampling/gaussian_combined.h>
#include <geosets/sampling/gaussian_rdhr.h>
#include <geosets/sampling/gaussian_rejection.h>
#include <geosets/sets/hpolytope.h>

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/vector.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace geosets::python::detail {

namespace nb = nanobind;
using DoubleArray = nb::ndarray<nb::numpy, double>;

template <typename SetType>
void validate_batch_shape(const std::vector<SetType*>& sets, const sampling::RowMajorMatrixXd& locs,
                          const sampling::RowMajorMatrixXd& scales) {
    // Only check invariants needed before indexing the batch. Each C++ sampler
    // owns validation of its set and Gaussian parameters.
    const Eigen::Index num_sets = static_cast<Eigen::Index>(sets.size());
    if (locs.rows() != num_sets || scales.rows() != num_sets) {
        throw std::invalid_argument(
            "Gaussian locs/scales row count must match the number of sets.");
    }
    for (const SetType* set : sets) {
        if (set == nullptr) {
            throw std::invalid_argument("Gaussian sampling received a null set pointer.");
        }
    }
}

inline DoubleArray pack_samples_to_numpy(const std::vector<sampling::RowMajorMatrixXd>& results,
                                         std::size_t samples_per_set, std::size_t dim) {
    const std::size_t num_sets = results.size();

    auto samples = std::make_unique<sampling::RowMajorMatrixXd>(
        static_cast<Eigen::Index>(samples_per_set * num_sets), static_cast<Eigen::Index>(dim));
    for (std::size_t sample_idx = 0; sample_idx < samples_per_set; ++sample_idx) {
        for (std::size_t set_idx = 0; set_idx < num_sets; ++set_idx) {
            const std::size_t output_row = sample_idx * num_sets + set_idx;
            samples->row(static_cast<Eigen::Index>(output_row)) =
                results[set_idx].row(static_cast<Eigen::Index>(sample_idx));
        }
    }

    double* data = samples->data();
    nb::capsule owner(samples.get(), [](void* pointer) noexcept {
        delete static_cast<sampling::RowMajorMatrixXd*>(pointer);
    });
    samples.release();
    return DoubleArray(data, {samples_per_set, num_sets, dim}, owner,
                       {static_cast<std::int64_t>(num_sets * dim), static_cast<std::int64_t>(dim),
                        static_cast<std::int64_t>(1)});
}

template <typename SetType, typename Sampler, typename... Args>
DoubleArray parallel_sample_sets(const std::vector<SetType*>& sets,
                                 const sampling::RowMajorMatrixXd& locs,
                                 const sampling::RowMajorMatrixXd& scales,
                                 std::size_t n_samples_per_set, Sampler sampler, Args... args) {
    validate_batch_shape(sets, locs, scales);
    const auto samples = geosets::parallel_execute(sets, [&](SetType* set, std::size_t set_idx) {
        const Eigen::Index row = static_cast<Eigen::Index>(set_idx);
        const Eigen::VectorXd loc = locs.row(row).transpose();
        const Eigen::VectorXd scale = scales.row(row).transpose();
        return sampler(*set, loc, scale, n_samples_per_set, args...);
    });
    return pack_samples_to_numpy(samples, n_samples_per_set, static_cast<std::size_t>(locs.cols()));
}

} // namespace geosets::python::detail

template <typename SetType>
void sampling_bindings(nanobind::module_& module, const std::string& typestr) {
    namespace nb = nanobind;

    // Single set

    module.def(("gaussian_rdhr_sampling_" + typestr).c_str(),
               geosets::sampling::gaussian_rdhr_samples, nb::arg("set"), nb::arg("loc"),
               nb::arg("scale"), nb::arg("n_samples") = 1, nb::arg("validate") = true);

    module.def(("gaussian_rejection_sampling_" + typestr).c_str(),
               geosets::sampling::gaussian_rejection_samples, nb::arg("set"), nb::arg("loc"),
               nb::arg("scale"), nb::arg("n_samples") = 1,
               nb::arg("rejection_limit") = geosets::sampling::DEFAULT_REJECTION_LIMIT,
               nb::arg("validate") = true);

    module.def(("gaussian_rejection_rdhr_fallback_sampling_" + typestr).c_str(),
               geosets::sampling::gaussian_rejection_rdhr_fallback_samples, nb::arg("set"),
               nb::arg("loc"), nb::arg("scale"), nb::arg("n_samples") = 1,
               nb::arg("rejection_limit") = 100, nb::arg("validate") = true);

    // Parallel for multiple set

    module.def(("batched_gaussian_rdhr_sampling_" + typestr).c_str(),
               [](const std::vector<SetType*>& sets,
                  const geosets::sampling::RowMajorMatrixXd& locs,
                  const geosets::sampling::RowMajorMatrixXd& scales, std::size_t n_samples_per_set,
                  bool validate) {
                   return geosets::python::detail::parallel_sample_sets(
                       sets, locs, scales, n_samples_per_set,
                       &geosets::sampling::gaussian_rdhr_samples, validate);
               },
               nb::arg("sets"), nb::arg("locs"), nb::arg("scales"), nb::arg("n_samples_per_set"),
               nb::arg("validate") = true);

    module.def(
        ("batched_gaussian_rejection_sampling_" + typestr).c_str(),
        [](const std::vector<SetType*>& sets, const geosets::sampling::RowMajorMatrixXd& locs,
           const geosets::sampling::RowMajorMatrixXd& scales, std::size_t n_samples_per_set,
           std::size_t rejection_limit, bool validate) {
            return geosets::python::detail::parallel_sample_sets(
                sets, locs, scales, n_samples_per_set,
                &geosets::sampling::gaussian_rejection_samples, rejection_limit, validate);
        },
        nb::arg("sets"), nb::arg("locs"), nb::arg("scales"), nb::arg("n_samples_per_set") = 1,
        nb::arg("rejection_limit") = geosets::sampling::DEFAULT_REJECTION_LIMIT,
        nb::arg("validate") = true);

    module.def(
        ("batched_gaussian_rejection_rdhr_fallback_sampling_" + typestr).c_str(),
        [](const std::vector<SetType*>& sets, const geosets::sampling::RowMajorMatrixXd& locs,
           const geosets::sampling::RowMajorMatrixXd& scales, std::size_t n_samples_per_set,
           std::size_t rejection_limit, bool validate) {
            return geosets::python::detail::parallel_sample_sets(
                sets, locs, scales, n_samples_per_set,
                &geosets::sampling::gaussian_rejection_rdhr_fallback_samples, rejection_limit,
                validate);
        },
        nb::arg("sets"), nb::arg("locs"), nb::arg("scales"), nb::arg("n_samples_per_set") = 1,
        nb::arg("rejection_limit") = 100, nb::arg("validate") = true);
}
