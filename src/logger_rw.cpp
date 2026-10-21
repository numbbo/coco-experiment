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
 * @brief The rw logger data type.
 *
 * @note Some fields from the observers (coco_observer as well as observer_rw) need to be copied here
 * because the observers can be deleted before the logger is finalized and we need these fields for
 * finalization.
 */
struct logger_rw_data_t {
  FILE *out_file;              /**< @brief File for logging. */
  size_t num_func_evaluations; /**< @brief The number of function evaluations performed so far. */
  size_t num_cons_evaluations; /**< @brief The number of evaluations of constraints performed so far. */

  double best_value;    /**< @brief The best-so-far value. */
  double current_value; /**< @brief The current value. */

  int log_vars;        /**< @brief Whether to log the decision values. */
  int log_cons;        /**< @brief Whether to log the constraints. */
  int log_only_better; /**< @brief Whether to log only solutions that are better than previous ones. */
  int log_time;        /**< @brief Whether to log evaluation time. */

  int precision_x;         /**< @brief Precision for outputting decision values. */
  int precision_f;         /**< @brief Precision for outputting objective values. */
  int precision_g;         /**< @brief Precision for outputting constraint values. */
  int log_discrete_as_int; /**< @brief Whether to output discrete variables in int or double format. */
};

static logger_rw_data_t *logger_rw_get_data(const coco_problem_t *problem);

/**
 * @brief Evaluates the function and constraints and outputs the information according to the
 * observer options.
 */
static void logger_rw_evaluate(coco_problem_t *problem, const double *x, double *y) {

  logger_rw_data_t *logger;
  coco_problem_t *inner_problem;
  std::vector<double> constraints(problem->number_of_constraints);
  size_t i;
  int log_this_time = 1;
  time_t start, end;

  logger = logger_rw_get_data(problem);
  inner_problem = coco_problem_transformed_get_inner_problem(problem);

  /* Time the evaluations */
  if (logger->log_time)
    time(&start);

  /* Evaluate the objective(s) */
  coco_evaluate_function(inner_problem, x, y);
  logger->num_func_evaluations++;

  if (problem->number_of_objectives == 1)
    logger->current_value = y[0];

  /* Evaluate the constraints */
  if (problem->number_of_constraints > 0)
    inner_problem->evaluate_constraint(x, constraints.data(), 0);
  logger->num_cons_evaluations = problem->evaluations_constraints;

  /* Time the evaluations */
  if (logger->log_time)
    time(&end);

  /* Log to the output file */
  if ((problem->number_of_objectives == 1) && (logger->current_value < logger->best_value))
    logger->best_value = logger->current_value;
  else if (problem->number_of_objectives == 1)
    log_this_time = !logger->log_only_better;
  if ((logger->num_func_evaluations == 1) || log_this_time) {
    fprintf(logger->out_file, "%lu\t", (unsigned long)logger->num_func_evaluations);
    fprintf(logger->out_file, "%lu\t", (unsigned long)logger->num_cons_evaluations);
    for (i = 0; i < problem->number_of_objectives; i++)
      fprintf(logger->out_file, "%+.*e\t", logger->precision_f, y[i]);
    if (logger->log_vars) {
      for (i = 0; i < problem->number_of_variables; i++) {
        if ((i < problem->number_of_integer_variables) && (logger->log_discrete_as_int))
          fprintf(logger->out_file, "%d\t", coco_double_to_int(x[i]));
        else
          fprintf(logger->out_file, "%+.*e\t", logger->precision_x, x[i]);
      }
    }
    if (logger->log_cons) {
      for (i = 0; i < problem->number_of_constraints; i++)
        fprintf(logger->out_file, "%+.*e\t", logger->precision_g, constraints[i]);
    }
    /* Log time in seconds */
    if (logger->log_time)
      fprintf(logger->out_file, "%.0f\t", difftime(end, start));
    fprintf(logger->out_file, "\n");
  }
  fflush(logger->out_file);
}

/**
 * @brief Frees the memory of the given rw logger.
 */
static void logger_rw_free(logger_rw_data_t *logger) {

  assert(logger != nullptr);
  if (logger->out_file != nullptr) {
    fclose(logger->out_file);
    logger->out_file = nullptr;
  }
}

/**
 * @brief The rw logger, which wraps the observed problem.
 */
struct logger_rw_t final : public coco_problem_transformed_t {
  logger_rw_data_t data;

  logger_rw_t(coco_problem_t *inner_problem, const char *observer_name)
      : coco_problem_transformed_t(inner_problem, observer_name), data() {}

  ~logger_rw_t() override { logger_rw_free(&data); }

  void evaluate_function(const double *x, double *y) override { logger_rw_evaluate(this, x, y); }
};

/**
 * @brief Returns the data of the logger.
 */
static logger_rw_data_t *logger_rw_get_data(const coco_problem_t *problem) {
  assert(dynamic_cast<const logger_rw_t *>(problem) != nullptr);
  return &static_cast<logger_rw_t *>(const_cast<coco_problem_t *>(problem))->data;
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
static coco_problem_t *logger_rw(coco_observer_t *observer, coco_problem_t *inner_problem) {

  logger_rw_t *problem;
  logger_rw_data_t *logger_data;
  observer_rw_t *observer_data;

  problem = new logger_rw_t(inner_problem, observer->observer_name.c_str());
  logger_data = &problem->data;
  logger_data->num_func_evaluations = 0;
  logger_data->num_cons_evaluations = 0;

  observer_data = observer_rw_get(observer);
  /* Copy values from the observes that you might need even if they do not exist any more */
  logger_data->precision_x = observer->precision_x;
  logger_data->precision_f = observer->precision_f;
  logger_data->precision_g = observer->precision_g;
  logger_data->log_discrete_as_int = observer->log_discrete_as_int;

  if (((observer_data->log_vars_mode == LOG_LOW_DIM) &&
       (inner_problem->number_of_variables > observer_data->low_dim_vars)) ||
      (observer_data->log_vars_mode == LOG_NEVER))
    logger_data->log_vars = 0;
  else
    logger_data->log_vars = 1;

  if (((observer_data->log_cons_mode == LOG_LOW_DIM) &&
       (inner_problem->number_of_constraints > observer_data->low_dim_cons)) ||
      (observer_data->log_cons_mode == LOG_NEVER) || (inner_problem->number_of_constraints == 0))
    logger_data->log_cons = 0;
  else
    logger_data->log_cons = 1;

  logger_data->log_only_better = (observer_data->log_only_better) && (inner_problem->number_of_objectives == 1);
  logger_data->log_time = observer_data->log_time;

  logger_data->best_value = DBL_MAX;
  logger_data->current_value = DBL_MAX;

  /* Construct file name */
  coco_create_directory(observer->result_folder);
  const std::string path_name =
      coco_join_path(observer->result_folder, coco_string_format("%s_rw.txt", coco_problem_get_id(inner_problem)));

  /* Open and initialize the output file */
  logger_data->out_file = fopen(path_name.c_str(), "a");
  if (logger_data->out_file == nullptr) {
    coco_error("logger_rw() failed to open file '%s'.", path_name.c_str());
    return nullptr; /* Never reached */
  }

  /* Output header information */
  fprintf(logger_data->out_file, "\n%% suite = '%s', problem_id = '%s', problem_name = '%s', coco_version = '%s'\n",
          coco_problem_get_suite(inner_problem)->suite_name.c_str(), coco_problem_get_id(inner_problem),
          coco_problem_get_name(inner_problem), coco_version);
  fprintf(logger_data->out_file, "%% f-evaluations | g-evaluations | %lu objective",
          (unsigned long)inner_problem->number_of_objectives);
  if (inner_problem->number_of_objectives > 1)
    fprintf(logger_data->out_file, "s");
  if (logger_data->log_vars)
    fprintf(logger_data->out_file, " | %lu variable", (unsigned long)inner_problem->number_of_variables);
  if (inner_problem->number_of_variables > 1)
    fprintf(logger_data->out_file, "s");
  if (logger_data->log_cons)
    fprintf(logger_data->out_file, " | %lu constraint", (unsigned long)inner_problem->number_of_constraints);
  if (inner_problem->number_of_constraints > 1)
    fprintf(logger_data->out_file, "s");
  if (logger_data->log_time)
    fprintf(logger_data->out_file, " | evaluation time (s)");
  fprintf(logger_data->out_file, "\n");


  return problem;
}
