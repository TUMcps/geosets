  # Constrained Gaussian Sampling in geosets

  ## 1. Sampling common infrastructure

  **Files**

  `include/geosets/sampling/common.h`

  `src/sampling/gaussian_common.cpp`


  - RowMajorMatrixXd: used for sampler outputs [samples_per_set, dim] and for packing the final C-contiguous Python output [samples_per_set, num_sets, dim].

  - validate_sample_inputs:

   if validate = False: only do structural checks: set dim > 0, loc size == dim, scale size == dim.

   if validate = True: additionally check loc/scale finite, scale positive, and set validity.

  - `sample_unit_uniform_open()` and `sample_standard_normal()` for random direction and 1d sample



  ## 2. RDHR sampler

  **Files**

  `include/geosets/sampling/gaussian_rdhr.h`

  `src/sampling/gaussian_rdhr.cpp`


  ### public API

  `gaussian_rdhr()` for one sample, output dim is [dim], and `gaussian_rdhr_samples()` for multiple samples per set, output dim is [samples_per_set, dim]

  ### detail helper

  `gaussian_rdhr_initial_point()` get the initial start point, either loc or center point, and `gaussian_rdhr_from_initial()` get the point sampled from rdhr.


  ### difference from origin acorl

 it has proposal_limit, which is a total attempt guard to avoid extreme long runtimes when RDHR proposals are repeatedly rejected by Metropolis filter or invalid.



  ## 3. Rejection sampler

  **Files**

  `include/geosets/sampling/gaussian_rejection.h`

  `src/sampling/gaussian_rejection.cpp`

  ### public API

 `gaussian_rejection_samples` returns [samples_per_set, dim], and throws if any sample is still not accepted after the initial draw and sampling_limit resampling rounds.

  ### internal helper
  `gaussian_rejection_samples_with_mask` return

   `struct RejectionResult {
      RowMajorMatrixXd samples;
      std::vector<std::uint8_t> not_accepted;
  };`

  the failure mask will be further used by RDHR fallback.

  `try_accept` will check if the point is inside the set



  ## 4. Rejection with RDHR fallback

  **Files**

  `include/geosets/sampling/gaussian_combined.h`

  `src/sampling/gaussian_combined.cpp`

  ### public API

  `gaussian_rejection_rdhr_fallback_samples()` do the final fallback using the rejection samples with mask from Rejection sampler and RDHR sampler. Only failed rejection samples are replaced by RDHR fallback; accepted rejection samples are kept unchanged.

  ## 5. HPolytope line intersection support

  **Files**

  - `include/geosets/sets/interface.h`
  - `include/geosets/sets/hpolytope.h`
  - `src/sets/hpolytope.cpp`

  `line_intersection_bounds` return the scalar bounds (lower, upper), which is used by RDHR

  `line_intersection` returns the boundary points {point + lower * direction, point + upper * direction}

  `line_intersection_bounds_impl` is the real function that calculate the bounds, which is now only implemented in HPolytope.


  ## 6. Python bindings and parallel execution

  **Files**

  - `bindings/python/sampling_bindings.h`
  - `bindings/python/nanobind_module.cpp`
  - `include/geosets/parallel.h`

  for `sampling_bindings.h`:

  1. validate_batch_shape which check if the size of locs and scales matches the sets.

  2. parallel_sample_sets uses geosets::parallel_execute. Each worker samples one set and returns [samples_per_set, dim]. pack_samples_to_numpy then converts these per-set outputs to [samples_per_set, num_sets, dim] for ACORL compatibility.

  3. In the end there are in total of 3 functions that are written by nanobind, which are RDHR, pure rejection and RDHR fallback

  for `parallel.h`:

  The reason of the change is that especially for pure rejection, only some of the workers will throw errors, which need to be catched separately.
