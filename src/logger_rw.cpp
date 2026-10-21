/**
 * @file logger_rw.cpp
 * @brief Implementation of the real-world logger.
 *
 * Can be used to log all (or just those that are better than the preceding) solutions with information
 * about objectives, decision variables (optional) and constraints (optional). See observer_rw() for
 * more information on the options. Produces one "txt" file for each problem function, dimension and
 * instance.
 *
 * @note This logger can be used with single- and multi-objective problems, but in the multi-objective
 * case, all solutions are always logged.
 */

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "coco.h"
#include "coco_internal.h"

#include "coco_utilities.cpp"
#include "coco_problem.cpp"
#include "coco_string.cpp"
#include "observer_rw.cpp"

namespace coco {

/**
 * @brief The rw logger, which wraps the observed problem.
 *
 * @note Some fields from the observers (coco_observer as well as observer_rw) are copied here because the
 * observers can be deleted before the logger is finalized and we need these fields for finalization.
 */
class LoggerRw final : public TransformedProblem {
public:
  LoggerRw(ObserverRw &observer, ProblemPtr inner);

  void EvaluateFunction(std::span<const double> x, std::span<double> y) override;

  /** @brief The decision values are not transformed. */
  bool Inverse(std::span<const double> inner_x, std::span<double> x) const override {
    return IdentityInverse(inner_x, x);
  }

private:
  File out_file_;                 /**< @brief File for logging. */
  size_t num_func_evaluations_{}; /**< @brief The number of function evaluations performed so far. */
  size_t num_cons_evaluations_{}; /**< @brief The number of evaluations of constraints performed so far. */
  double best_found_value_{};     /**< @brief The best-so-far value. */
  double current_value_{};        /**< @brief The current value. */
  bool log_vars_{};               /**< @brief Whether to log the decision values. */
  bool log_cons_{};               /**< @brief Whether to log the constraints. */
  bool log_only_better_{};        /**< @brief Whether to log only solutions that are better than previous ones. */
  bool log_time_{};               /**< @brief Whether to log evaluation time. */
  int precision_x_{};             /**< @brief Precision for outputting decision values. */
  int precision_f_{};             /**< @brief Precision for outputting objective values. */
  int precision_g_{};             /**< @brief Precision for outputting constraint values. */
  bool log_discrete_as_int_{};    /**< @brief Whether to output discrete variables in int or double format. */
};

/**
 * @brief Evaluates the function and constraints and outputs the information according to the
 * observer options.
 */
void LoggerRw::EvaluateFunction(std::span<const double> x, std::span<double> y) {
  std::vector<double> constraints(number_of_constraints_);
  bool log_this_time = true;
  time_t start, end;

  /* Time the evaluations */
  if (log_time_)
    time(&start);

  /* Evaluate the objective(s) */
  inner_problem_->EvaluateFunction(x, y);
  num_func_evaluations_++;

  if (number_of_objectives_ == 1)
    current_value_ = y[0];

  /* Evaluate the constraints */
  if (number_of_constraints_ > 0)
    inner_problem_->EvaluateConstraint(x, constraints, 0);
  num_cons_evaluations_ = evaluations_constraints();

  /* Time the evaluations */
  if (log_time_)
    time(&end);

  /* Log to the output file */
  if ((number_of_objectives_ == 1) && (current_value_ < best_found_value_))
    best_found_value_ = current_value_;
  else if (number_of_objectives_ == 1)
    log_this_time = !log_only_better_;
  if ((num_func_evaluations_ == 1) || log_this_time) {
    fprintf(out_file_, "%zu\t", num_func_evaluations_);
    fprintf(out_file_, "%zu\t", num_cons_evaluations_);
    for (size_t i = 0; i < number_of_objectives_; i++)
      fprintf(out_file_, "%+.*e\t", precision_f_, y[i]);
    if (log_vars_) {
      for (size_t i = 0; i < number_of_variables_; i++) {
        if ((i < number_of_integer_variables_) && (log_discrete_as_int_))
          fprintf(out_file_, "%d\t", DoubleToInt(x[i]));
        else
          fprintf(out_file_, "%+.*e\t", precision_x_, x[i]);
      }
    }
    if (log_cons_) {
      for (size_t i = 0; i < number_of_constraints_; i++)
        fprintf(out_file_, "%+.*e\t", precision_g_, constraints[i]);
    }
    /* Log time in seconds */
    if (log_time_)
      fprintf(out_file_, "%.0f\t", difftime(end, start));
    fprintf(out_file_, "\n");
  }
  fflush(out_file_);
}

/**
 * @brief Initializes the rw logger.
 *
 * Copies all observer field values that are needed after initialization into logger field values for two
 * reasons:
 * - If the observer is deleted before the suite, the observer is not available anymore when the logger
 * is finalized.
 * - This reduces function calls.
 */
LoggerRw::LoggerRw(ObserverRw &observer, ProblemPtr inner)
    : TransformedProblem(std::move(inner), observer.observer_name) {
  num_func_evaluations_ = 0;
  num_cons_evaluations_ = 0;
  /* Copy values from the observes that you might need even if they do not exist any more */
  precision_x_ = observer.precision_x;
  precision_f_ = observer.precision_f;
  precision_g_ = observer.precision_g;
  log_discrete_as_int_ = observer.log_discrete_as_int;

  if (((observer.log_vars_mode_ == LOG_LOW_DIM) && (inner_problem_->number_of_variables() > observer.low_dim_vars_)) ||
      (observer.log_vars_mode_ == LOG_NEVER))
    log_vars_ = 0;
  else
    log_vars_ = 1;

  if (((observer.log_cons_mode_ == LOG_LOW_DIM) &&
       (inner_problem_->number_of_constraints() > observer.low_dim_cons_)) ||
      (observer.log_cons_mode_ == LOG_NEVER) || (inner_problem_->number_of_constraints() == 0))
    log_cons_ = 0;
  else
    log_cons_ = 1;

  log_only_better_ = (observer.log_only_better_) && (inner_problem_->number_of_objectives() == 1);
  log_time_ = observer.log_time_;

  best_found_value_ = DBL_MAX;
  current_value_ = DBL_MAX;

  /* Construct file name */
  CreateDirectory(observer.result_folder);
  const std::string path_name =
      JoinPath(observer.result_folder, StringFormat("%s_rw.txt", coco_problem_get_id(inner_problem_.get())));

  /* Open and initialize the output file */
  out_file_.Open(path_name, "a");
  if (out_file_ == nullptr) {
    coco_error("LoggerRw() failed to open file '%s'.", path_name.c_str());
    return; /* Never reached */
  }

  /* Output header information */
  fprintf(out_file_, "\n%% suite = '%s', problem_id = '%s', problem_name = '%s', coco_version = '%s'\n",
          inner_problem_->suite()->name().c_str(), coco_problem_get_id(inner_problem_.get()),
          coco_problem_get_name(inner_problem_.get()), coco_version);
  fprintf(out_file_, "%% f-evaluations | g-evaluations | %zu objective", inner_problem_->number_of_objectives());
  if (inner_problem_->number_of_objectives() > 1)
    fprintf(out_file_, "s");
  if (log_vars_)
    fprintf(out_file_, " | %zu variable", inner_problem_->number_of_variables());
  if (inner_problem_->number_of_variables() > 1)
    fprintf(out_file_, "s");
  if (log_cons_)
    fprintf(out_file_, " | %zu constraint", inner_problem_->number_of_constraints());
  if (inner_problem_->number_of_constraints() > 1)
    fprintf(out_file_, "s");
  if (log_time_)
    fprintf(out_file_, " | evaluation time (s)");
  fprintf(out_file_, "\n");
}

ProblemPtr ObserverRw::AddLogger(ProblemPtr problem) { return std::make_unique<LoggerRw>(*this, std::move(problem)); }

} // namespace coco
