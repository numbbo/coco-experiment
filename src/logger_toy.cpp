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

/**
 * @brief The toy logger, which wraps the observed problem.
 */
class logger_toy_t final : public coco_problem_transformed_t {
public:
  logger_toy_t(observer_toy_t &observer, coco_problem_ptr_t inner);

  void evaluate_function(std::span<const double> x, std::span<double> y) override;

  /** @brief The decision values are not transformed. */
  bool inverse(std::span<const double> inner_x, std::span<double> x) const override {
    return identity_inverse(inner_x, x);
  }

private:
  FILE *log_file{}; /**< @brief Pointer to the file already prepared for logging (owned by the observer). */
  coco_observer_log_targets_t targets; /**< @brief Triggers based on logarithmic target values. */
  size_t number_of_evaluations{};      /**< @brief The number of evaluations performed so far. */
  int precision_x{};                   /**< @brief Precision for outputting decision values. */
  int precision_f{};                   /**< @brief Precision for outputting objective values. */
};

/**
 * @brief Evaluates the function, increases the number of evaluations and outputs information based on the
 * targets that have been hit.
 */
void logger_toy_t::evaluate_function(std::span<const double> x, std::span<double> y) {

  inner_problem->evaluate_function(x, y);
  number_of_evaluations++;

  /* Output the solution when a new target that has been hit */
  if (targets.trigger(y[0])) {
    fprintf(log_file, "%zu\t%.*e\t%.*e", number_of_evaluations, precision_f, y[0], precision_f, targets.get_value());
    for (size_t i = 0; i < number_of_variables; i++) {
      fprintf(log_file, "\t%.*e", precision_x, x[i]);
    }
    fprintf(log_file, "\n");
  }

  /* Flush output so that impatient users can see the progress */
  fflush(log_file);
}

/**
 * @brief Initializes the toy logger.
 */
logger_toy_t::logger_toy_t(observer_toy_t &observer, coco_problem_ptr_t inner)
    : coco_problem_transformed_t(std::move(inner), observer.observer_name) {
  if (inner_problem->get_number_of_objectives() != 1) {
    coco_warning("logger_toy_t(): The toy logger shouldn't be used to log a problem with %zu objectives",
                 inner_problem->get_number_of_objectives());
  }

  /* Initialize the logger */
  number_of_evaluations = 0;
  targets = coco_observer_log_targets_t(observer.number_target_triggers, observer.log_target_precision);
  log_file = observer.log_file;
  precision_x = observer.precision_x;
  precision_f = observer.precision_f;

  /* Output initial information */
  assert(inner_problem->get_suite());
  fprintf(log_file, "\n");
  fprintf(log_file, "suite = '%s', problem_id = '%s', problem_name = '%s', coco_version = '%s'\n",
          inner_problem->get_suite()->get_name().c_str(), coco_problem_get_id(inner_problem.get()),
          coco_problem_get_name(inner_problem.get()), coco_version);
  fprintf(log_file, "%% evaluation number | function value | target hit | %zu variables \n",
          inner_problem->get_number_of_variables());
}

coco_problem_ptr_t observer_toy_t::add_logger(coco_problem_ptr_t problem) {
  return std::make_unique<logger_toy_t>(*this, std::move(problem));
}
