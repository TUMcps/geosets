#include <cstdio>
#include <geosets/parallel.h>
#include <geosets/sets/hpolytope.h>
#include <geosets/sets/vpolytope.h>

// #ifdef GS_USE_OPENMP
// #include <omp.h>
// #define OMP_PARALLEL_FOR _Pragma("omp parallel for")
// #else
// #define OMP_PARALLEL_FOR
// #endif

int main() {
    OMP_PARALLEL_FOR
    for (int i = 0; i < 10; ++i) {
        std::printf("Hello from iteration %d\n", i);

#ifdef GS_USE_OPENMP
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        std::printf("  (OpenMP enabled: thread %d of %d)\n", tid, nthreads);
#endif

        auto a = geosets::HPolytope::from_unit_box(2);
        // std::printf("%s\n\n", a.str().c_str());
    }

    std::vector<geosets::HPolytope> polys;
    for (int i = 0; i < 10; ++i)
        polys.emplace_back(geosets::HPolytope::from_unit_box(3));

    geosets::set_omp_num_threads(8);
    auto centers =
        parallel_execute(polys, [](const geosets::HPolytope& p, size_t) { return p.center(); });
    std::printf("%zu", centers.size());
}
