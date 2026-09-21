#pragma once

#include <atomic>
#include <exception>
#include <random>
#include <type_traits>
#include <vector>

#ifdef GS_USE_OPENMP
#include <omp.h>
#define OMP_PARALLEL_FOR _Pragma("omp parallel for")
#else
#define OMP_PARALLEL_FOR
#endif

namespace geosets {

namespace detail {
inline std::atomic<unsigned int> g_rng_seed{5489u};   // mt19937 default seed
inline std::atomic<unsigned int> g_rng_generation{0}; // bumped on each reseed
} // namespace detail

// Per-thread RNG: each thread mutates its own engine, so parallel_execute workers draw without
// contention or data races (std::rand()/a shared engine would either lock or race per draw).
inline std::mt19937& thread_rng() {
    static std::atomic<unsigned int> next_index{0};
    thread_local unsigned int index = next_index.fetch_add(1, std::memory_order_relaxed);
    thread_local unsigned int seen_gen = ~0u;
    thread_local std::mt19937 engine;
    const unsigned int gen = detail::g_rng_generation.load(std::memory_order_relaxed);
    if (seen_gen != gen) {
        std::seed_seq seq{detail::g_rng_seed.load(std::memory_order_relaxed), index};
        engine.seed(seq);
        seen_gen = gen;
    }
    return engine;
}

// Reseed the RNG base. Live per-thread engines pick up the new seed on their next draw.
inline void seed_rng(unsigned int seed) {
    detail::g_rng_seed.store(seed, std::memory_order_relaxed);
    detail::g_rng_generation.fetch_add(1, std::memory_order_relaxed);
}

// Function to define the number of threads used
inline void set_omp_num_threads(int nthreads) {
#ifdef GS_USE_OPENMP
    omp_set_num_threads(nthreads);
#else
    (void)nthreads;
    std::printf("OpenMP not enabled; set_omp_num_threads has no effect.\n");
#endif
}

// Generic parallel execution for returning results
template <typename SetType, typename Func>
auto parallel_execute(const std::vector<SetType>& sets, Func func)
    -> std::vector<decltype(func(std::declval<SetType>(), size_t{}))> {
    using ResultType = decltype(func(std::declval<SetType>(), size_t{}));
    // std::vector<bool> packs bits, so concurrent writes to neighbouring elements race.
    using StoredType = std::conditional_t<std::is_same_v<ResultType, bool>, char, ResultType>;
    std::vector<StoredType> results(sets.size());
    std::vector<std::exception_ptr> errors(sets.size());

    OMP_PARALLEL_FOR
    for (size_t i = 0; i < sets.size(); ++i) {
        try {
            results[i] = func(sets[i], i);
        } catch (...) {
            errors[i] = std::current_exception();
        }
    }

    for (const auto& error : errors) {
        if (error) std::rethrow_exception(error);
    }
    if constexpr (std::is_same_v<StoredType, ResultType>) {
        return results;
    } else {
        return std::vector<ResultType>(results.begin(), results.end());
    }
}

// Generic parallel execution for functions with void return type
template <typename SetType, typename Func>
void parallel_execute_void(std::vector<SetType>& sets, Func func) {

    // Debug
    // #ifdef GS_USE_OPENMP
    //     std::printf("Running parallel_execute with %d threads.\n", omp_get_max_threads());
    // #else
    //     std::printf("Running parallel_execute without OpenMP.\n");
    // #endif

    OMP_PARALLEL_FOR
    for (size_t i = 0; i < sets.size(); ++i) {
        func(sets[i], i);
    }
}

} // namespace geosets
