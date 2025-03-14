/**
 * @file logger_bbob.c
 * @brief Implementation of the bbob logger.
 *
 * Logs the performance of an optimizer on single-objective problems with or without constraints
 * and with or without knowing their true optimal value.
 *
 * It produces a number of files:
 * - The "info" files contain high-level information on the performed experiment. One .info file
 * is created for each function and contains information on all the problem instances for that
 * function.
 * - The remaining files ("dat", "tdat", "rdat" and "mdat" files) contain detailed information
 * on the performance of an optimizer. Logging in the "dat" files is triggered when performance
 * targets are reached (they are either uniform in logarithmic or linear scale, depending on
 * whether the problem has a known optimal value or not*1, respectively), logging in the "tdat"
 * files is triggered by the number of function evaluations, logging in the "rdat" files is
 * triggered by the restart of an algorithm and logging in the "mdat" is triggered every time
 * a solution is recommended and the performance targets are reached (the same targets are 
 * used as for the "dat" files).
 *
 * *1 Note that the decision about which performance targets to use is done on the suite level
 * and passed to the logger during initialization.
 */

#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <errno.h>

#include "coco.h"

#include "coco_utilities.c"
#include "coco_problem.c"
#include "coco_string.c"
#include "coco_observer.c"
#include "observer_bbob.c"

/**
 * @brief Factor used in logged indicator (f-f*)^+ + sum_i g_i^+ in front of the sum
 */
#define LOGGER_BBOB_WEIGHT_CONSTRAINTS 1e0

/**
 * Penalty used when logging infeasible solutions of problems with unknown optimal value
 */
#define LOGGER_BBOB_INFEASIBLE_PENALTY 1e10

/**
 * @brief The bbob logger data type.
 */
typedef struct {
  coco_observer_t *observer; /**< @brief Pointer to the observer (might be NULL at the end) */
  char *suite_name;          /**< @brief The suite name */
  int is_initialized;        /**< @brief Whether the logger was already initialized */
  int algorithm_restarted;   /**< @brief Whether the algorithm has restarted (output information to .rdat file). */

  FILE *info_file; /**< @brief Index file */
  FILE *dat_file;  /**< @brief File with function value aligned data */
  FILE *tdat_file; /**< @brief File with number of evaluations aligned data */
  FILE *rdat_file; /**< @brief File with restart information */
  FILE *mdat_file; /**< @brief File with evaluated recommendations */

  size_t num_func_evaluations;        /**< @brief The number of function evaluations performed so far. */
  size_t num_cons_evaluations;        /**< @brief The number of evaluations of constraints performed so far. */
  size_t last_logged_eval_dat;        /**< @brief The last evaluation logged in the dat file (needed for finalization) */
  size_t last_logged_eval_tdat;       /**< @brief The last evaluation logged in the tdat file (needed for finalization) */
  size_t last_logged_eval_mdat;       /**< @brief The last evaluation logged in the mdat file (needed for finalization) */

  double *best_found_solution;        /**< @brief The best found solution to this problem. */
  double best_found_value;            /**< @brief The best found value for this problem. */
  double current_value;               /**< @brief The current value for this problem. */
  double optimal_value;               /**< @brief The optimal value for this problem (if it is known, otherwise a reference value). */

  double *best_recommended_solution;  /**< @brief The best recommended solution for this problem. */
  double best_recommended_value;      /**< @brief The best recommended value for this problem. */

  size_t function;                    /**< @brief Suite-dependent function number */
  size_t instance;                    /**< @brief Suite-dependent instance number */
  size_t number_of_variables;         /**< @brief Number of all variables */
  size_t number_of_integer_variables; /**< @brief Number of integer variables */
  size_t number_of_constraints;       /**< @brief Number of constraints */
  int log_discrete_as_int;            /**< @brief Whether to output discrete variables in integer or double format. */

  coco_observer_targets_t *targets;         /**< @brief Triggers based on target values. */
  coco_observer_evaluations_t *evaluations; /**< @brief Triggers based on the number of evaluations. */
  coco_observer_targets_t *mdat_targets;    /**< @brief Triggers based on target values for recommendations. */

} logger_bbob_data_t;

/**
 * @brief The data format used by the logger
 *
 * Back to 5 columns, 5-th column writes single digit constraint values.
 *
 * Previous formats:
 *
 * -> "old":
 * function evaluation |
 * noise-free fitness - Fopt (7.948000000000e+01) |
 * best noise-free fitness - Fopt |
 * measured fitness |
 * best measured fitness |
 * x1 | x2...
 *
 * -> "bbob":
 * f evaluations |
 * g evaluations |
 * best noise-free fitness - Fopt |
 * noise-free fitness - Fopt (%13.12e) |
 * measured fitness |
 * best measured fitness |
 * x1 | x2...
 */
static const char *logger_bbob_data_format = "bbob-new2";

/**
 * @brief The header used by the logger in .?dat files
 *
 * See also logger_bbob_data_format
 */
static const char *logger_bbob_header = "%% "
                                        "f evaluations | "
                                        "g evaluations | "
                                        "best noise-free fitness - %s (%13.12e) + sum g_i+ | "
                                        "measured fitness | "
                                        "best measured fitness or single-digit g-values | "
                                        "x1 | "
                                        "x2...\n";

/**
 * @brief Discretized constraint value, ~8 + log10(c), in a single digit.
 *
 * -\infty..0 -> 0
 *    0..1e-7 -> 1
 * 1e-7..1e-6 -> 2
 *    ...
 * 1e-1..1    -> 8
 *   >1       -> 9
 */
static int logger_bbob_single_digit_constraint(const double c) {
  const double limits[9] = {0, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1};
  int i;

  for (i = 0; i < 9; ++i)
    if (c <= limits[i])
      return i;
  return 9;
}

/**
 * @brief Decides whether a new line in the info file (and therefore a new set of .?dat files)
 * should be started
 *
 * A new line should be started when the current dimension or the current function differ from
 * the last logged dimension and function.
 */
static int logger_bbob_start_new_line(coco_observer_t *observer, size_t current_dimension, size_t current_function) {
  observer_bbob_data_t *observer_data = ((observer_bbob_data_t *)observer->data);
  assert(observer_data != NULL);

  if ((current_function == observer_data->last_function) && (current_dimension == observer_data->last_dimension))
    return 0;
  return 1;
}

/**
 * @brief Outputs a formatted line to a data file
 */
static void logger_bbob_output(FILE *data_file, logger_bbob_data_t *logger, const double *x, double current_value,
                               double best_value, const double *constraints) {
  /* This function contains many hard-coded values (10.9, 22, 5.4) that could be read through
   * observer options */
  size_t i;

  fprintf(data_file, "%lu %lu %+10.9e %+10.9e ", (unsigned long)logger->num_func_evaluations,
          (unsigned long)logger->num_cons_evaluations, best_value - logger->optimal_value, current_value);

  if ((logger->number_of_constraints > 0) && (constraints != NULL)) {
    for (i = 0; i < logger->number_of_constraints; ++i) {
      /* print 01234567890123..., may happen in the last line of .tdat */
      fprintf(data_file, "%d", constraints ? logger_bbob_single_digit_constraint(constraints[i]) : (int)(i % 10));
    }
  } else {
    fprintf(data_file, "%+10.9e", best_value);
  }

  if (logger->number_of_variables < 22) {
    for (i = 0; i < logger->number_of_variables; i++) {
      if ((i < logger->number_of_integer_variables) && (logger->log_discrete_as_int))
        fprintf(data_file, " %d", coco_double_to_int(x[i]));
      else
        fprintf(data_file, " %+5.4e", x[i]);
    }
  }
  fprintf(data_file, "\n");

  /* Flush output so that impatient users can see progress.
   * Otherwise it can take a long time until the output appears.
   */
  fflush(data_file);
}

/**
 * @brief Opens the file in append mode
 */
static void logger_bbob_open_file(FILE **file, const char *file_path) {
  if (*file == NULL) {
    *file = fopen(file_path, "a");
    if (*file == NULL) {
      coco_error("logger_bbob_open_file(): Error opening file: %s\nError: %d", file_path, errno);
    }
  }
}

/**
 * @brief Creates the data file (if it didn't exist before) and opens it
 */
static void logger_bbob_open_data_file(FILE **data_file, const char *path, const char *file_name,
                                       const char *file_extension) {
  char file_path[COCO_PATH_MAX + 2] = {0};
  char relative_file_path[COCO_PATH_MAX + 2] = {0};
  strncpy(relative_file_path, file_name, COCO_PATH_MAX - strlen(relative_file_path) - 1);
  strncat(relative_file_path, file_extension, COCO_PATH_MAX - strlen(relative_file_path) - 1);
  coco_join_path(file_path, sizeof(file_path), path, relative_file_path, NULL);
  logger_bbob_open_file(data_file, file_path);
}

/**
 * @brief Creates the info file (if it didn't exist before) and opens it
 */
static void logger_bbob_open_info_file(logger_bbob_data_t *logger, const char *folder, const char *function_string,
                                       const char *data_file_name, const char *suite_name, int start_new_line) {
  char data_file_path[COCO_PATH_MAX + 2] = {0};
  int add_empty_line = 0;
  char file_name[COCO_PATH_MAX + 2] = {0};
  char file_path[COCO_PATH_MAX + 2] = {0};
  FILE **info_file;
  FILE *tmp_file;
  coco_observer_t *observer;
  observer_bbob_data_t *observer_data;

  coco_debug("Started logger_bbob_open_info_file()");

  assert(logger != NULL);
  assert(logger->observer != NULL);
  observer = (coco_observer_t *)(logger->observer);
  observer_data = (observer_bbob_data_t *)(observer->data);
  assert(observer_data != NULL);

  strncpy(data_file_path, data_file_name, COCO_PATH_MAX - strlen(data_file_path) - 1);

  info_file = &(logger->info_file);

  strncpy(file_name, observer_data->prefix, COCO_PATH_MAX - strlen(file_name) - 1);
  strncat(file_name, "_f", COCO_PATH_MAX - strlen(file_name) - 1);
  strncat(file_name, function_string, COCO_PATH_MAX - strlen(file_name) - 1);
  strncat(file_name, ".info", COCO_PATH_MAX - strlen(file_name) - 1);
  coco_join_path(file_path, sizeof(file_path), folder, file_name, NULL);

  if (*info_file == NULL) {
    add_empty_line = 0;
    /* If the file already exists, an empty line is needed */
    tmp_file = fopen(file_path, "r");
    if (tmp_file != NULL) {
      add_empty_line = 1;
      fclose(tmp_file);
    }
    logger_bbob_open_file(info_file, file_path);
    if (start_new_line) {
      if (add_empty_line)
        fprintf(*info_file, "\n");

      fprintf(*info_file,
              "suite = '%s', funcId = %lu, DIM = %lu, Precision = %.3e, "
              "algId = '%s', coco_version = '%s', logger = '%s', "
              "data_format = '%s', settings = '%s'\n",
              suite_name, (unsigned long)logger->function, (unsigned long)logger->number_of_variables, pow(10, -8),
              observer->algorithm_name, coco_version, observer->observer_name,
              logger_bbob_data_format, observer->settings);
      fprintf(*info_file, "%% %s\n", observer->algorithm_info);
      /* data_file_path does not have the extension */
      fprintf(*info_file, "%s.dat", data_file_path);
    }
  }
  coco_debug("Ended   logger_bbob_open_info_file()");
}

/**
 * @brief Generates the files and folders needed by the logger if they don't already exist
 */
static void logger_bbob_initialize(logger_bbob_data_t *logger, int is_opt_known) {

  char relative_path[COCO_PATH_MAX + 2] = {0}; /* Relative path to the .dat file from where the .info file is */
  char folder_path[COCO_PATH_MAX + 2] = {0};
  char *function_string;
  char *dimension_string;
  char *str_pointer;
  char *relative_path_pointer = NULL;
  int start_new_line;
  observer_bbob_data_t *observer_data;

  coco_debug("Started logger_bbob_initialize()");

  if (is_opt_known == 0)
    str_pointer = "Fref";
  else
    str_pointer = "Fopt";

  assert(logger != NULL);
  assert(logger->observer != NULL);
  observer_data = ((observer_bbob_data_t *)((coco_observer_t *)logger->observer)->data);

  function_string = coco_strdupf("%lu", (unsigned long)logger->function);
  dimension_string = coco_strdupf("%lu", (unsigned long)logger->number_of_variables);

  /* Prepare paths and names */
  strncpy(relative_path, "data_f", COCO_PATH_MAX);
  strncat(relative_path, function_string, COCO_PATH_MAX - strlen(relative_path) - 1);
  coco_join_path(folder_path, sizeof(folder_path), logger->observer->result_folder, relative_path, NULL);
  coco_create_directory(folder_path);
  strncat(relative_path, "/bbobexp_f", COCO_PATH_MAX - strlen(relative_path) - 1);
  strncat(relative_path, function_string, COCO_PATH_MAX - strlen(relative_path) - 1);
  strncat(relative_path, "_DIM", COCO_PATH_MAX - strlen(relative_path) - 1);
  strncat(relative_path, dimension_string, COCO_PATH_MAX - strlen(relative_path) - 1);

  /* Find out whether a new line in the info file (and a new .?dat file) is needed */
  start_new_line = logger_bbob_start_new_line(logger->observer, logger->number_of_variables, logger->function);
  if (start_new_line) {
    relative_path_pointer = coco_strdup(relative_path);
    coco_create_unique_filename(logger->observer->result_folder, &relative_path_pointer, ".dat");
    strncpy(observer_data->last_dat_file, relative_path_pointer, COCO_PATH_MAX);
    observer_data->last_function = logger->function;
    observer_data->last_dimension = logger->number_of_variables;
  } else {
    relative_path_pointer = coco_strdup(observer_data->last_dat_file);
  }

  /* info file */
  logger_bbob_open_info_file(logger, logger->observer->result_folder, function_string, relative_path_pointer,
                             logger->suite_name, start_new_line);
  fprintf(logger->info_file, ", %lu", (unsigned long)logger->instance);

  /* data files */
  logger_bbob_open_data_file(&(logger->dat_file), logger->observer->result_folder, relative_path_pointer, ".dat");
  fprintf(logger->dat_file, logger_bbob_header, str_pointer, logger->optimal_value);
  logger_bbob_open_data_file(&(logger->tdat_file), logger->observer->result_folder, relative_path_pointer, ".tdat");
  fprintf(logger->tdat_file, logger_bbob_header, str_pointer, logger->optimal_value);
  logger_bbob_open_data_file(&(logger->rdat_file), logger->observer->result_folder, relative_path_pointer, ".rdat");
  fprintf(logger->rdat_file, logger_bbob_header, str_pointer, logger->optimal_value);
  logger_bbob_open_data_file(&(logger->mdat_file), logger->observer->result_folder, relative_path_pointer, ".mdat");
  fprintf(logger->mdat_file, logger_bbob_header, str_pointer, logger->optimal_value);

  logger->is_initialized = 1;
  coco_free_memory(dimension_string);
  coco_free_memory(function_string);
  coco_free_memory(relative_path_pointer);

  coco_debug("Ended   logger_bbob_initialize()");
}

/**
 * @brief Evaluates the function, increases the number of evaluations and outputs information according to
 * observer options.
 */
static void logger_bbob_evaluate(coco_problem_t *problem, const double *x, double *y) {
  size_t i;
  double y_logged, max_value = 0, sum_constraints;
  double *constraints = NULL;
  logger_bbob_data_t *logger = (logger_bbob_data_t *)coco_problem_transformed_get_data(problem);
  coco_problem_t *inner_problem = coco_problem_transformed_get_inner_problem(problem);
  const int is_feasible = problem->number_of_constraints <= 0 || coco_is_feasible(inner_problem, x, NULL);

  if (problem->is_tainted || inner_problem->is_tainted) {
    coco_error("Though shalt not cheat! Trying to evaluate a tainted observed problem. Do _not_ call "
               "coco_problem_best_parameter() on a problem you want to benchmark.");
  }

  coco_debug("Started logger_bbob_evaluate()");

  if (!logger->is_initialized) {
    logger_bbob_initialize(logger, problem->is_opt_known);
  }
  if ((coco_log_level >= COCO_DEBUG) && logger->num_func_evaluations == 0) {
    coco_debug("%4lu: ", (unsigned long)inner_problem->suite_dep_index);
    coco_debug("on problem %s ... ", coco_problem_get_id(inner_problem));
  }

  /* Fulfill contract of a COCO evaluate function */
  coco_evaluate_function(inner_problem, x, y);
  logger->num_func_evaluations++;
  logger->current_value = y[0];
  y_logged = y[0];
  if (inner_problem->is_noisy == 1) {
    logger->current_value = inner_problem->last_noise_free_values[0];
    y_logged = inner_problem->last_noise_free_values[0];
  }
  if (coco_is_nan(y_logged))
    y_logged = NAN_FOR_LOGGING;
  else if (coco_is_inf(y_logged))
    y_logged = INFINITY_FOR_LOGGING;

  /* Do sanity check */
  if ((problem->is_opt_known) && (is_feasible)) {
    /* Infeasible solutions can have much better y0 values */
    assert(y_logged + 1e-13 >= logger->optimal_value);
  }

  /* Evaluate the constraints */
  if (problem->number_of_constraints > 0) {
    constraints = coco_allocate_vector(problem->number_of_constraints);
    inner_problem->evaluate_constraint(inner_problem, x, constraints, 0);
  }
  logger->num_cons_evaluations = problem->evaluations_constraints;

  /* Compute the sum of positive constraint values */
  sum_constraints = 0;
  for (i = 0; i < problem->number_of_constraints; ++i) {
    if (constraints[i] > 0)
      sum_constraints += constraints[i];
  }
  sum_constraints *= LOGGER_BBOB_WEIGHT_CONSTRAINTS; /* do this before the checks */
  if (coco_is_nan(sum_constraints))
    sum_constraints = NAN_FOR_LOGGING;
  else if (coco_is_inf(sum_constraints))
    sum_constraints = INFINITY_FOR_LOGGING;

  if (problem->is_opt_known)
    max_value = coco_double_max(y_logged, logger->optimal_value);
  else {
    /* Problems with unknown optimal values */
    if (sum_constraints > 0) {
      max_value = LOGGER_BBOB_INFEASIBLE_PENALTY; /* Infeasible solution */
    } else
      max_value = y_logged; /* Feasible solution */
  }

  /* Update logger state
   * At logger->number_of_evaluations == 1 the logger->best_found_value is not initialized,
   * also compare to max_value to not potentially be thrown off by weird values in y[0]
   */
  if (logger->num_func_evaluations == 1 || (max_value + sum_constraints < logger->best_found_value)) {
    logger->best_found_value = max_value + sum_constraints;
    for (i = 0; i < problem->number_of_variables; i++)
      logger->best_found_solution[i] = x[i]; /* May well be infeasible */

    /* Add a line in the .dat file is the performance target is reached by a feasible solution
     * and always at the first evaluation */
    if (logger->num_func_evaluations == 1 ||
        coco_observer_targets_trigger(logger->targets, logger->best_found_value - logger->optimal_value)) {
      logger_bbob_output(logger->dat_file, logger, x, y_logged, logger->best_found_value, constraints);
      logger->last_logged_eval_dat = logger->num_func_evaluations;
    }
  }

  /* Add a line in the .tdat file for specific number of evaluations*/
  if (coco_observer_evaluations_trigger(logger->evaluations,
                                        logger->num_func_evaluations + logger->num_cons_evaluations)) {
    logger_bbob_output(logger->tdat_file, logger, x, y_logged, logger->best_found_value, constraints);
    logger->last_logged_eval_tdat = logger->num_func_evaluations;
  }

  /* Add a line in the .rdat file if the algorithm was restarted */
  if (logger->algorithm_restarted) {
    logger_bbob_output(logger->rdat_file, logger, x, y_logged, logger->best_found_value, constraints);
    logger->algorithm_restarted = 0;
  }

  /* Free allocated memory */
  if (problem->number_of_constraints > 0)
    coco_free_memory(constraints);

  coco_debug("Ended   logger_bbob_evaluate()");
}

/**
 * @brief Evaluates the function and outputs information according to observer options to the file with
 * recommendations. The evaluation result is not returned and the evaluation counter is not increased.
 */
static void logger_bbob_recommend(coco_problem_t *problem, const double *x) {
  double y_logged;
  double *constraints = NULL, *y = NULL;
  logger_bbob_data_t *logger = (logger_bbob_data_t *)coco_problem_transformed_get_data(problem);
  coco_problem_t *inner_problem = coco_problem_transformed_get_inner_problem(problem);
  const int is_feasible = problem->number_of_constraints <= 0 || coco_is_feasible(inner_problem, x, NULL);

  if (!logger->is_initialized) {
    logger_bbob_initialize(logger, problem->is_opt_known);
  }

  /* Fulfill contract of a COCO evaluate function, but do not increase the evaluation counters */
  y = coco_allocate_vector(problem->number_of_objectives);
  coco_evaluate_function(inner_problem, x, y);

  y_logged = y[0];
  if (coco_is_nan(y_logged))
    y_logged = NAN_FOR_LOGGING;
  else if (coco_is_inf(y_logged))
    y_logged = INFINITY_FOR_LOGGING;
  coco_free_memory(y);

  /* Do sanity check */
  if ((problem->is_opt_known) && (is_feasible)) {
    /* Infeasible solutions can have much better y0 values */
    assert(y_logged + 1e-13 >= logger->optimal_value);
  }

  /* Evaluate the constraints */
  if (problem->number_of_constraints > 0) {
    constraints = coco_allocate_vector(problem->number_of_constraints);
    inner_problem->evaluate_constraint(inner_problem, x, constraints, 0);
  }

  /* Compute the sum of positive constraint values */
  double sum_constraints;
  sum_constraints = 0.0;
  for (size_t i = 0; i < problem->number_of_constraints; ++i) {
    if (constraints[i] > 0)
      sum_constraints += constraints[i];
  }
  sum_constraints *= LOGGER_BBOB_WEIGHT_CONSTRAINTS; /* do this before the checks */
  if (coco_is_nan(sum_constraints))
    sum_constraints = NAN_FOR_LOGGING;
  else if (coco_is_inf(sum_constraints))
    sum_constraints = INFINITY_FOR_LOGGING;

  double max_value;
  if (problem->is_opt_known) 
    max_value = coco_double_max(y_logged, logger->optimal_value);
  else {
    /* Problems with unknown optimal values */
    if (sum_constraints > 0) {
      max_value = LOGGER_BBOB_INFEASIBLE_PENALTY; /* Infeasible solution */
    } else
      max_value = y_logged; /* Feasible solution */
  }

  /* Update recommender state
   * The logger->best_recommended_value is initialized to INFINITY_FOR_LOGGING.
   */
  if (max_value + sum_constraints < logger->best_recommended_value) {
    for (size_t i = 0; i < problem->number_of_variables; i++)
      logger->best_recommended_solution[i] = x[i]; /* May well be infeasible */

    logger->best_recommended_value = max_value + sum_constraints;
    /* Add a line in the .mdat file is the performance target is reached by a feasible solution */
    if (coco_observer_targets_trigger(logger->mdat_targets, logger->best_recommended_value - logger->optimal_value)) {
      logger_bbob_output(logger->mdat_file, logger, x, y_logged, logger->best_recommended_value, constraints);
      logger->last_logged_eval_mdat = logger->num_func_evaluations;
    }
  }

  /* Free allocated memory */
  if (problem->number_of_constraints > 0)
    coco_free_memory(constraints);
}

/**
 * @brief Takes care of some final output and frees the logger data
 */
static void logger_bbob_free(void *stuff) {
  assert(stuff != NULL);

  logger_bbob_data_t *logger = (logger_bbob_data_t *)stuff;
  coco_observer_t *observer; /* The observer might not exist at this point */

  coco_debug("Started logger_bbob_free()");

  if ((coco_log_level >= COCO_DEBUG) && logger && logger->num_func_evaluations > 0) {
    coco_debug("best f=%e after %lu fevals (done observing)\n", logger->best_found_value,
               (unsigned long)logger->num_func_evaluations);
  }
  if (logger->info_file != NULL) {
    fprintf(logger->info_file, ":%lu|%.1e", (unsigned long)logger->num_func_evaluations,
            logger->best_found_value - logger->optimal_value);
    fclose(logger->info_file);
    logger->info_file = NULL;
  }

  if (logger->dat_file != NULL) {
    if (logger->num_func_evaluations != logger->last_logged_eval_dat)
      logger_bbob_output(logger->dat_file, logger, logger->best_found_solution, logger->best_found_value, 
                         logger->best_found_value, NULL);
    fclose(logger->dat_file);
    logger->dat_file = NULL;
  }

  if (logger->tdat_file != NULL) {
    if (logger->num_func_evaluations != logger->last_logged_eval_tdat)
      logger_bbob_output(logger->tdat_file, logger, logger->best_found_solution, logger->best_found_value, 
                         logger->best_found_value, NULL);
    fclose(logger->tdat_file);
    logger->tdat_file = NULL;
  }

  if (logger->rdat_file != NULL) {
    fclose(logger->rdat_file);
    logger->rdat_file = NULL;
  }

  if (logger->mdat_file != NULL) {
    if ((logger->num_func_evaluations != logger->last_logged_eval_mdat) && (logger->last_logged_eval_mdat != 0))
      logger_bbob_output(logger->mdat_file, logger, logger->best_recommended_solution, logger->best_recommended_value, 
                         logger->best_recommended_value, NULL);
    fclose(logger->mdat_file);
    logger->mdat_file = NULL;
  }

  if (logger->best_found_solution != NULL) {
    coco_free_memory(logger->best_found_solution);
    logger->best_found_solution = NULL;
  }

  if (logger->best_recommended_solution != NULL) {
    coco_free_memory(logger->best_recommended_solution);
    logger->best_recommended_solution = NULL;
  }

  if (logger->targets != NULL) {
    coco_observer_targets_free(logger->targets);
    logger->targets = NULL;
  }

  if (logger->evaluations != NULL) {
    coco_observer_evaluations_free(logger->evaluations);
    logger->evaluations = NULL;
  }

  if (logger->mdat_targets != NULL) {
    coco_observer_targets_free(logger->mdat_targets);
    logger->mdat_targets = NULL;
  }

  observer = logger->observer;
  if ((observer != NULL) && (observer->is_active == 1)) {
    if (observer->data != NULL) {
      /* If the observer still exists (if it does not, observed_problem does not matter any longer) */
      ((observer_bbob_data_t *)observer->data)->observed_problem = NULL;
    }
  }

  coco_debug("Ended   logger_bbob_free()");
}

/*
 * @brief Sets the pointer to the observer to NULL
 *
 * Has to be taken care of here due to
 */
static void logger_bbob_data_nullify_observer(void *stuff) {
  logger_bbob_data_t *logger = (logger_bbob_data_t *)stuff;
  logger->observer = NULL;
}

/**
 * @brief Saves the information that the algorithm has restarted
 */
static void logger_bbob_signal_restart(coco_problem_t *problem) {

  logger_bbob_data_t *logger = (logger_bbob_data_t *)coco_problem_transformed_get_data(problem);
  assert(logger);

  if (logger->num_func_evaluations > 0)
    logger->algorithm_restarted = 1;
}

static coco_problem_t *logger_bbob(coco_observer_t *observer, coco_problem_t *inner_problem) {
  coco_problem_t *problem;
  logger_bbob_data_t *logger_data;
  observer_bbob_data_t *observer_data;
  coco_suite_t *suite;
  size_t i;

  coco_debug("Started logger_bbob()");

  assert(inner_problem);
  assert(inner_problem->suite);
  suite = (coco_suite_t *)inner_problem->suite;
  if (inner_problem->number_of_objectives != 1) {
    coco_warning("logger_bbob(): The bbob logger shouldn't be used to log a problem with %d objectives",
                 inner_problem->number_of_objectives);
  }

  assert(observer != NULL);
  observer_data = (observer_bbob_data_t *)observer->data;
  assert(observer_data != NULL);
  if (observer_data->observed_problem != NULL) {
    coco_error("logger_bbob(): The observed problem must be closed before a new problem can be observed");
  }

  logger_data = (logger_bbob_data_t *)coco_allocate_memory(sizeof(*logger_data));
  logger_data->observer = observer;
  logger_data->suite_name = coco_problem_get_suite(inner_problem)->suite_name;
  logger_data->is_initialized = 0;
  logger_data->algorithm_restarted = 0;

  logger_data->info_file = NULL;
  logger_data->dat_file = NULL;
  logger_data->tdat_file = NULL;
  logger_data->rdat_file = NULL;
  logger_data->mdat_file = NULL;

  logger_data->num_func_evaluations = 0;
  logger_data->num_cons_evaluations = 0;
  logger_data->last_logged_eval_mdat = 0;
  /* logger_data->last_logged_evaluation = 0; */

  logger_data->best_found_solution = coco_allocate_vector(inner_problem->number_of_variables);
  logger_data->best_found_value = DBL_MAX;
  logger_data->optimal_value = *(inner_problem->best_value);

  logger_data->current_value = DBL_MAX;
  
  logger_data->best_recommended_solution = coco_allocate_vector(inner_problem->number_of_variables);
  logger_data->best_recommended_value = INFINITY_FOR_LOGGING;

  logger_data->function = coco_problem_get_suite_dep_function(inner_problem);
  logger_data->instance = coco_problem_get_suite_dep_instance(inner_problem);
  logger_data->number_of_variables = inner_problem->number_of_variables;
  logger_data->number_of_integer_variables = inner_problem->number_of_integer_variables;
  logger_data->log_discrete_as_int = observer->log_discrete_as_int;
  logger_data->number_of_constraints = inner_problem->number_of_constraints;

  /* Initialize triggers performance and evaluation triggers depending on the observer option */
  logger_data->targets = coco_observer_targets(suite->known_optima, observer->lin_target_precision,
                                               observer->number_target_triggers, observer->log_target_precision);
  logger_data->evaluations =
      coco_observer_evaluations(observer->base_evaluation_triggers, inner_problem->number_of_variables);
  logger_data->mdat_targets = coco_observer_targets(suite->known_optima, observer->lin_target_precision,
                                                    observer->number_target_triggers, observer->log_target_precision);

  coco_debug("Ended   logger_bbob()");

  problem = coco_problem_transformed_allocate(inner_problem, logger_data, logger_bbob_free, observer->observer_name);
  problem->evaluate_function = logger_bbob_evaluate;
  problem->recommend_solution = logger_bbob_recommend;

  observer_data->observed_problem = problem;
  return problem;
}
