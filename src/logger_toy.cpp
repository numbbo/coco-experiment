/**
 * @file logger_toy.cpp
 * @brief Implementation of the toy logger.
 *
 * Logs the evaluation number, function value the target hit and all the variables each time a target has
 * been hit.
 */

#include <stdio.h>
#include <assert.h>

#include "coco.h"

#include "coco_utilities.cpp"
#include "coco_problem.cpp"
#include "coco_string.cpp"
#include "observer_toy.cpp"

namespace coco {

/**
 * @brief The toy logger, which wraps the observed problem.
 */
class LoggerToy final : public Logger {
public:
  LoggerToy(ObserverToy &observer, ProblemPtr inner);

  void EvaluateFunction(std::span<const double> x, std::span<double> y) override;

  /** @brief The decision values are not transformed. */
  bool Inverse(std::span<const double> inner_x, std::span<double> x) const override {
    return IdentityInverse(inner_x, x);
  }

private:
  FILE *log_file_{}; /**< @brief Pointer to the file already prepared for logging (owned by the observer). */
  ObserverLogTargets targets_;          /**< @brief Triggers based on logarithmic target values. */
  size_t number_of_evaluations_{};      /**< @brief The number of evaluations performed so far. */
  int precision_x_{};                   /**< @brief Precision for outputting decision values. */
  int precision_f_{};                   /**< @brief Precision for outputting objective values. */
};

/**
 * @brief Evaluates the function, increases the number of evaluations and outputs information based on the
 * targets that have been hit.
 */
void LoggerToy::EvaluateFunction(std::span<const double> x, std::span<double> y) {

  inner_problem_->EvaluateFunction(x, y);
  number_of_evaluations_++;

  /* Output the solution when a new target that has been hit */
  if (targets_.Trigger(y[0])) {
    fprintf(log_file_, "%zu\t%.*e\t%.*e", number_of_evaluations_, precision_f_, y[0], precision_f_, targets_.value());
    for (size_t i = 0; i < number_of_variables_; i++) {
      fprintf(log_file_, "\t%.*e", precision_x_, x[i]);
    }
    fprintf(log_file_, "\n");
  }

  /* Flush output so that impatient users can see the progress */
  fflush(log_file_);
}

/**
 * @brief Initializes the toy logger.
 */
LoggerToy::LoggerToy(ObserverToy &observer, ProblemPtr inner) : Logger(observer, std::move(inner)) {
  if (inner_problem_->number_of_objectives() != 1) {
    coco_warning("LoggerToy(): The toy logger shouldn't be used to log a problem with %zu objectives",
                 inner_problem_->number_of_objectives());
  }

  /* Initialize the logger */
  number_of_evaluations_ = 0;
  targets_ = ObserverLogTargets(observer.number_target_triggers, observer.log_target_precision);
  log_file_ = observer.log_file_;
  precision_x_ = observer.precision_x;
  precision_f_ = observer.precision_f;

  /* Output initial information */
  assert(inner_problem_->suite());
  fprintf(log_file_, "\n");
  fprintf(log_file_, "suite = '%s', problem_id = '%s', problem_name = '%s', coco_version = '%s'\n",
          inner_problem_->suite()->name().c_str(), coco_problem_get_id(inner_problem_.get()),
          coco_problem_get_name(inner_problem_.get()), coco_version);
  fprintf(log_file_, "%% evaluation number | function value | target hit | %zu variables \n",
          inner_problem_->number_of_variables());
}

ProblemPtr ObserverToy::AddLogger(ProblemPtr problem) { return std::make_unique<LoggerToy>(*this, std::move(problem)); }

} // namespace coco
