#pragma once

#include <stdexcept>
#include <string>

namespace geosets {

struct NotImplemented : public std::logic_error {
  public:
    NotImplemented() : std::logic_error("Feature not yet implemented") {}
    NotImplemented(std::string message) : std::logic_error(message) {}
};

struct InfiniteSetException : std::runtime_error {
    InfiniteSetException() : std::runtime_error("Operation not defined for infinite sets") {}
    InfiniteSetException(std::string message) : std::runtime_error(message) {}
};

struct EmptySetException : std::runtime_error {
    EmptySetException() : std::runtime_error("Operation not defined for empty sets") {}
    EmptySetException(std::string message) : std::runtime_error(message) {}
};

struct SetOperationException : std::runtime_error {
    using std::runtime_error::runtime_error;
};

struct DimensionMismatchException : SetOperationException {
    DimensionMismatchException(size_t expected, size_t got)
        : SetOperationException("Dimension mismatch: expected " + std::to_string(expected) +
                                ", got " + std::to_string(got)) {}
};

struct FailedOptimizationProblem : SetOperationException {
    FailedOptimizationProblem()
        : SetOperationException("The optimization problem failed to solve.") {}
    FailedOptimizationProblem(std::string message) : SetOperationException(std::move(message)) {}
};
struct InfeasibleOptimizationProblem : SetOperationException {
    InfeasibleOptimizationProblem()
        : SetOperationException("The optimization problem is infeasible.") {}
};
struct NotOptimalOptimizationProblem : SetOperationException {
    NotOptimalOptimizationProblem()
        : SetOperationException("The optimization problem could not find an optimal solution.") {}
};

struct DegenerateSetException : SetOperationException {
    DegenerateSetException() : SetOperationException("The set is degenerate.") {}
};

} // namespace geosets
