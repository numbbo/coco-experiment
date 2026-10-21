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

/**
 * @brief The rw logger, which wraps the observed problem.
 *
 * @note Some fields from the observers (coco_observer as well as observer_rw) are copied here because the
 * observers can be deleted before the logger is finalized and we need these fields for finalization.
 */
class logger_rw_t final : public coco_problem_transformed_t {
public:
  logger_rw_t(observer_rw_t &observer, coco_problem_ptr_t inner);

  void evaluate_function(std::span<const double> x, std::span<double> y) override;

  /** @brief The decision values are not transformed. */
  bool inverse(std::span<const double> inner_x, std::span<double> x) const override {
    return identity_inverse(inner_x, x);
  }

private:
  coco_file_t out_file;          /**< @brief File for logging. */
  size_t num_func_evaluations{}; /**< @brief The number of function evaluations performed so far. */
  size_t num_cons_evaluations{}; /**< @brief The number of evaluations of constraints performed so far. */
  double best_found_value{};     /**< @brief The best-so-far value. */
  double current_value{};        /**< @brief The current value. */
  bool log_vars{};               /**< @brief Whether to log the decision values. */
  bool log_cons{};               /**< @brief Whether to log the constraints. */
  bool log_only_better{};        /**< @brief Whether to log only solutions that are better than previous ones. */
  bool log_time{};               /**< @brief Whether to log evaluation time. */
  int precision_x{};             /**< @brief Precision for outputting decision values. */
  int precision_f{};             /**< @brief Precision for outputting objective values. */
  int precision_g{};             /**< @brief Precision for outputting constraint values. */
  bool log_discrete_as_int{};    /**< @brief Whether to output discrete variables in int or double format. */
};

/**
 * @brief Evaluates the function and constraints and outputs the information according to the
 * observer options.
 */
void logger_rw_t::evaluate_function(std::span<const double> x, std::span<double> y) {
  std::vector<double> constraints(number_of_constraints);
  bool log_this_time = true;
  time_t start, end;

  /* Time the evaluations */
  if (log_time)
    time(&start);

  /* Evaluate the objective(s) */
  inner_problem->evaluate_function(x, y);
  num_func_evaluations++;

  if (number_of_objectives == 1)
    current_value = y[0];

  /* Evaluate the constraints */
  if (number_of_constraints > 0)
    inner_problem->evaluate_constraint(x, constraints, 0);
  num_cons_evaluations = get_evaluations_constraints();

  /* Time the evaluations */
  if (log_time)
    time(&end);

  /* Log to the output file */
  if ((number_of_objectives == 1) && (current_value < best_found_value))
    best_found_value = current_value;
  else if (number_of_objectives == 1)
    log_this_time = !log_only_better;
  if ((num_func_evaluations == 1) || log_this_time) {
    fprintf(out_file, "%zu\t", num_func_evaluations);
    fprintf(out_file, "%zu\t", num_cons_evaluations);
    for (size_t i = 0; i < number_of_objectives; i++)
      fprintf(out_file, "%+.*e\t", precision_f, y[i]);
    if (log_vars) {
      for (size_t i = 0; i < number_of_variables; i++) {
        if ((i < number_of_integer_variables) && (log_discrete_as_int))
          fprintf(out_file, "%d\t", coco_double_to_int(x[i]));
        else
          fprintf(out_file, "%+.*e\t", precision_x, x[i]);
      }
    }
    if (log_cons) {
      for (size_t i = 0; i < number_of_constraints; i++)
        fprintf(out_file, "%+.*e\t", precision_g, constraints[i]);
    }
    /* Log time in seconds */
    if (log_time)
      fprintf(out_file, "%.0f\t", difftime(end, start));
    fprintf(out_file, "\n");
  }
  fflush(out_file);
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
logger_rw_t::logger_rw_t(observer_rw_t &observer, coco_problem_ptr_t inner)
    : coco_problem_transformed_t(std::move(inner), observer.observer_name) {
  num_func_evaluations = 0;
  num_cons_evaluations = 0;
  /* Copy values from the observes that you might need even if they do not exist any more */
  precision_x = observer.precision_x;
  precision_f = observer.precision_f;
  precision_g = observer.precision_g;
  log_discrete_as_int = observer.log_discrete_as_int;

  if (((observer.log_vars_mode == LOG_LOW_DIM) && (inner_problem->get_number_of_variables() > observer.low_dim_vars)) ||
      (observer.log_vars_mode == LOG_NEVER))
    log_vars = 0;
  else
    log_vars = 1;

  if (((observer.log_cons_mode == LOG_LOW_DIM) &&
       (inner_problem->get_number_of_constraints() > observer.low_dim_cons)) ||
      (observer.log_cons_mode == LOG_NEVER) || (inner_problem->get_number_of_constraints() == 0))
    log_cons = 0;
  else
    log_cons = 1;

  log_only_better = (observer.log_only_better) && (inner_problem->get_number_of_objectives() == 1);
  log_time = observer.log_time;

  best_found_value = DBL_MAX;
  current_value = DBL_MAX;

  /* Construct file name */
  coco_create_directory(observer.result_folder);
  const std::string path_name =
      coco_join_path(observer.result_folder, coco_string_format("%s_rw.txt", coco_problem_get_id(inner_problem.get())));

  /* Open and initialize the output file */
  out_file.open(path_name, "a");
  if (out_file == nullptr) {
    coco_error("logger_rw_t() failed to open file '%s'.", path_name.c_str());
    return; /* Never reached */
  }

  /* Output header information */
  fprintf(out_file, "\n%% suite = '%s', problem_id = '%s', problem_name = '%s', coco_version = '%s'\n",
          inner_problem->get_suite()->get_name().c_str(), coco_problem_get_id(inner_problem.get()),
          coco_problem_get_name(inner_problem.get()), coco_version);
  fprintf(out_file, "%% f-evaluations | g-evaluations | %zu objective", inner_problem->get_number_of_objectives());
  if (inner_problem->get_number_of_objectives() > 1)
    fprintf(out_file, "s");
  if (log_vars)
    fprintf(out_file, " | %zu variable", inner_problem->get_number_of_variables());
  if (inner_problem->get_number_of_variables() > 1)
    fprintf(out_file, "s");
  if (log_cons)
    fprintf(out_file, " | %zu constraint", inner_problem->get_number_of_constraints());
  if (inner_problem->get_number_of_constraints() > 1)
    fprintf(out_file, "s");
  if (log_time)
    fprintf(out_file, " | evaluation time (s)");
  fprintf(out_file, "\n");
}

coco_problem_ptr_t observer_rw_t::add_logger(coco_problem_ptr_t problem) {
  return std::make_unique<logger_rw_t>(*this, std::move(problem));
}
