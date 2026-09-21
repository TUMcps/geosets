#include <geosets/exceptions.h>
#include <geosets/highs_wrapper.h>

namespace geosets::highs_wrapper {

HighsModelStatus solve_lp(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                          const Eigen::VectorXd& objective, Eigen::VectorXd* solution,
                          bool check_scaling) {
    const int m = static_cast<int>(A.rows());
    const int n = static_cast<int>(A.cols());

    // Thread-local instance: avoids per-call construction overhead. HiGHS is
    // not internally thread-safe, so each thread keeps its own. Highs is
    // non-copyable/non-movable, so default-construct in place and init once.
    thread_local Highs highs;
    thread_local const bool highs_initialized = [] {
        highs.setOptionValue("output_flag", false);
        highs.setOptionValue("log_to_console", false);
        highs.setOptionValue("presolve", "off"); // Faster for the small LP problems here
        return true;
    }();
    (void)highs_initialized;
    highs.clearModel();

    // HiGHS ignores matrix coefficients with |a_ij| <= small_matrix_value (mirrored by the fill
    // loop below). Read it rather than hardcode, so we track whatever HiGHS actually applies.
    double small_matrix_value = 0.0;
    highs.getOptionValue("small_matrix_value", small_matrix_value);

    // Opt-in guard: a constraint row entirely at/below that threshold would be silently culled and
    // the LP solved without it, giving a wrong or failed result on a well-defined but ill-scaled
    // set. Reject it with a clear message instead of a cryptic HiGHS "Solve error".
    if (check_scaling) {
        for (int i = 0; i < m; ++i) {
            const double row_max = A.row(i).cwiseAbs().maxCoeff();
            if (row_max > 0.0 && row_max <= small_matrix_value) {
                throw FailedOptimizationProblem(
                    "LP constraint row " + std::to_string(i) +
                    " has all coefficients at or below HiGHS's small_matrix_value tolerance; "
                    "the set is too ill-scaled to solve reliably.");
            }
        }
    }

    HighsLp lp;
    lp.num_col_ = n;
    lp.num_row_ = m;

    lp.col_lower_.assign(n, -kHighsInf);
    lp.col_upper_.assign(n, kHighsInf);

    lp.col_cost_.resize(n);
    for (int i = 0; i < n; ++i) {
        lp.col_cost_[i] = -objective(i);
    }
    lp.offset_ = 0.0;

    lp.row_lower_.assign(m, -kHighsInf);
    lp.row_upper_.resize(m);
    for (int i = 0; i < m; ++i) {
        lp.row_upper_[i] = b(i);
    }

    HighsInt nnz = 0;
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.start_.resize(static_cast<size_t>(n + 1));
    for (int j = 0; j < n; ++j) {
        lp.a_matrix_.start_[static_cast<size_t>(j)] = nnz;
        for (int i = 0; i < m; ++i) {
            const double val = A(i, j);
            // Don't include very small values (they would give a warning in HighsStatus)
            if (std::abs(val) > small_matrix_value) {
                lp.a_matrix_.index_.push_back(i);
                lp.a_matrix_.value_.push_back(val);
                ++nnz;
            }
        }
    }
    lp.a_matrix_.start_[static_cast<size_t>(n)] = nnz;

    const HighsStatus pass_status = highs.passModel(lp);
    if (pass_status != HighsStatus::kOk) {
        throw FailedOptimizationProblem();
    }

    const HighsStatus run_status = highs.run();
    if (run_status != HighsStatus::kOk) {
        throw FailedOptimizationProblem();
    }

    const HighsModelStatus model_status = highs.getModelStatus();
    if (model_status == HighsModelStatus::kOptimal && solution != nullptr) {
        const HighsSolution& sol = highs.getSolution();
        solution->resize(n);
        for (int i = 0; i < n; ++i) {
            (*solution)(i) = sol.col_value[static_cast<size_t>(i)];
        }
    }

    return model_status;
}

} // namespace geosets::highs_wrapper
