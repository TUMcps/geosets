# The CPS geometric sets library

![pipeline](https://gitlab.lrz.de/cps/geosets/badges/main/pipeline.svg)
![coverage](https://gitlab.lrz.de/cps/geosets/badges/main/coverage.svg)

## Installation

### Dependencies
No C/C++ dependencies are required, since they are fetched directly when not installed on the system. However, for **faster compile times**, please install the following dependencies
- `highs`
- `boost`
- `eigen`
- `qhull`

On mac:
```sh
brew install highs boost eigen qhull
```

### Python uv installation
Add the gitlab repo to your uv sources, and then add `geosets-py` as a dependency.

```toml
dependencies = [
    "geosets-py"
]
[tool.uv.sources]
geosets_py = { git = "ssh://git@gitlab.lrz.de/cps/geosets.git", branch = "main" }
```

### Python pip installation
```sh
pip install git+ssh://git@gitlab.lrz.de/cps/geosets.git@main
```

**Important note on omp:**
Other python package might ship a bundled `omp` library (e.g., `pytorch`). Enabling `omp` for `geosets` then, results in two different `omp` libraries at runtime, which causes undefined behavior and crashes.
The cmake logic in `omp.cmake` can account for this by trying to find the bundled `omp` library in the python environment. Consequently, for this to work, `geosets-py` needs to be compiled from source inside the python venv, after the other package which bundles `omp`. The only way to ensure this is to exclude `geoesets-py` from the dependencies, and add it later manually with
```sh
export CMAKE_ARGS="-DGS_USE_OPENMP=On" pip install geosets-py
```



### C++ CMake installation
Simply include this in your cmake scripts

```cmake
include(FetchContent)
FetchContent_Declare(
    geosets
    GIT_REPOSITORY git@gitlab.lrz.de:cps/geosets.git
    GIT_TAG        main
)
FetchContent_MakeAvailable(geosets)
```

And then link geosets against your target
```cmake
target_link_libraries(your_target_name geosets::geosets)
```


### Dev installation
The project uses `uv`
```sh
uv sync
```

## Implemented operations
::include{file=documentation/table.md}

## Constrained Gaussian sampling

Constrained Gaussian sampling currently supports `HPolytope`. The Python bindings expose
direct RDHR, pure rejection, and rejection with position-level RDHR fallback:

```python
samples = geosets_py.batched_gaussian_rejection_rdhr_fallback_sampling_hpolytope(
    sets, locs, scales, samples_per_set, sampling_limit
)
```

`locs` and `scales` have shape `[num_sets, dim]`. The batched Python API defaults to
`validate=False`; use `validate=True` for untrusted sets or distribution parameters.
Structural checks needed to index the batch and sampler safely remain enabled.
All three samplers return a NumPy array with shape
`[samples_per_set, num_sets, dim]`. Sampling runs in parallel across sets; the binding
assembles the sample-major ACORL-compatible layout only when creating the NumPy result.


### Development
The Weirdest bug ever. Make sure to use `.eval()` at the end when creating random vectors with Eigen.
```C++
Eigen::VectorXd::Random(dim).eval();
```

## Parallelization
If you want to use parallelized execution with `OpenMP`, set the CMake flag `-DGS_USE_OPENMP=ON` when configuring the build.
For the python bindings you can enable the flag in the `pyproject.toml` file.
There are general functions in `parallel.h` that can be used for parallelized execution of operations, for example
```C++
std::vector<geosets::HPolytope> polys;
for (int i = 0; i < 10; ++i)
    polys.emplace_back(geosets::HPolytope::from_unit_box(3));

geosets::set_omp_num_threads(8);
auto centers = parallel_execute(polys, [](const geosets::HPolytope& p, size_t) { return p.center(); });
```

### Performance
Testing showed that using the batched execution in C++ only provides marginal speedup compared to using loops in python. The real improvement only comes when enabling `OpenMP`!
```python
# Python loop
for i in range(batch_size):
    center = sets[i].center(validate=False)
    boundary_point = sets[i].boundary_point(points[i], validate=False)

# C++ batch
centers = gsp.batched_center_hpolytope(sets)
boundary_points = gsp.batched_boundary_point_hpolytope(sets, points)

# Results for dim=3, 128 sets without omp enabled
# Python computation time: 0.0072 seconds
# C++ computation time: 0.0061 seconds
```

### cddlib thread safety
If you want to use parallelization for functions that use cddlib under the hood, make sure to have at least version 0.94n installed. Previous versions are not thread safe (https://github.com/cddlib/cddlib/releases).
