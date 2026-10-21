/**
 * @file logger_bbob.cpp
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

#include "coco_utilities.cpp"
#include "coco_problem.cpp"
#include "coco_string.cpp"
#include "coco_observer.cpp"
#include "observer_bbob.cpp"

/**
 * @brief Factor used in logged indicator (f-f*)^+ + sum_i g_i^+ in front of the sum
 */
#define LOGGER_BBOB_WEIGHT_CONSTRAINTS 1e0

/**
 * Penalty used when logging infeasible solutions of problems with unknown optimal value
 */
#define LOGGER_BBOB_INFEASIBLE_PENALTY 1e10

/**
 * Record variable values only when #variables is below this threshold, was: 22
 */
#define LOGGER_BBOB_NUMBER_OF_VARIABLES_THRESHOLD 7

/**
 * @brief The bbob logger, which wraps the observed problem.
 */
struct logger_bbob_t final : public coco_problem_transformed_t {
  observer_bbob_t *observer{}; /**< @brief Pointer to the observer (might be NULL at the end) */
  std::string suite_name;    /**< @brief The suite name */
  int is_initialized{};      /**< @brief Whether the logger was already initialized */
  int algorithm_restarted{}; /**< @brief Whether the algorithm has restarted (output information to .rdat file). */

  FILE *info_file{}; /**< @brief Index file */
  FILE *dat_file{};  /**< @brief File with function value aligned data */
  FILE *tdat_file{}; /**< @brief File with number of evaluations aligned data */
  FILE *rdat_file{}; /**< @brief File with restart information */
  FILE *mdat_file{}; /**< @brief File with evaluated recommendations */

  size_t num_func_evaluations{};  /**< @brief The number of function evaluations performed so far. */
  size_t num_cons_evaluations{};  /**< @brief The number of evaluations of constraints performed so far. */
  size_t last_logged_eval_dat{};  /**< @brief The last evaluation logged in the dat file (needed for finalization) */
  size_t last_logged_eval_tdat{}; /**< @brief The last evaluation logged in the tdat file (needed for finalization) */
  size_t last_logged_eval_mdat{}; /**< @brief The last evaluation logged in the mdat file (needed for finalization) */

  std::vector<double> best_found_solution; /**< @brief The best found solution to this problem. */
  double best_found_value{};               /**< @brief The best found value for this problem. */
  double current_value{};                  /**< @brief The current value for this problem. */
  double
      optimal_value{}; /**< @brief The optimal value for this problem (if it is known, otherwise a reference value). */

  std::vector<double> best_recommended_solution; /**< @brief The best recommended solution for this problem. */
  double best_recommended_value{};               /**< @brief The best recommended value for this problem. */

  int log_discrete_as_int{}; /**< @brief Whether to output discrete variables in integer or double format. */

  coco_observer_targets_t targets;      /**< @brief Triggers based on target values. */
  coco_observer_evaluations_t evaluation_triggers; /**< @brief Triggers based on the number of evaluations. */
  coco_observer_targets_t mdat_targets; /**< @brief Triggers based on target values for recommendations. */

  logger_bbob_t(observer_bbob_t &observer, coco_problem_ptr_t inner);

  ~logger_bbob_t() override;

  void evaluate_function(const double *x, double *y) override;
  void recommend_solution(const double *x) override;

  /** @brief Saves the information that the algorithm has restarted. */
  void signal_restart();

  /** @brief Returns the bbob logger behind the given problem (which must be one). */
  static logger_bbob_t *from(const coco_problem_t *problem);

private:
  /** @brief Generates the files and folders needed by the logger if they don't already exist. */
  void initialize();

  /** @brief Creates the info file (if it didn't exist before) and opens it. */
  void open_info_file(const std::string &folder, const std::string &function_string, const std::string &data_file_name,
                      int start_new_line);

  /** @brief Outputs a formatted line to a data file. */
  void output(FILE *data_file, const double *x, double logged_value, double logged_best_value,
              const double *constraints) const;
};

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
static int logger_bbob_start_new_line(const observer_bbob_t *observer, size_t current_dimension,
                                      size_t current_function) {
  assert(observer != nullptr);

  if ((current_function == observer->last_function) && (current_dimension == observer->last_dimension))
    return 0;
  return 1;
}

/**
 * @brief Outputs a formatted line to a data file
 */
void logger_bbob_t::output(FILE *data_file, const double *x, double logged_value, double logged_best_value,
                           const double *constraints) const {
  /* This function contains many hard-coded values (10.9, 22, 5.4) that could be read through
   * observer options */
  size_t i;

  fprintf(data_file, "%lu %lu %+10.9e %+10.9e ", (unsigned long)num_func_evaluations,
          (unsigned long)num_cons_evaluations, logged_best_value - optimal_value, logged_value);

  if ((number_of_constraints > 0) && (constraints != nullptr)) {
    for (i = 0; i < number_of_constraints; ++i) {
      /* print 01234567890123..., may happen in the last line of .tdat */
      fprintf(data_file, "%d", constraints ? logger_bbob_single_digit_constraint(constraints[i]) : (int)(i % 10));
    }
  } else {
    fprintf(data_file, "%+10.9e", logged_best_value);
  }

  if (number_of_variables < LOGGER_BBOB_NUMBER_OF_VARIABLES_THRESHOLD) {
    for (i = 0; i < number_of_variables; i++) {
      if ((i < number_of_integer_variables) && (log_discrete_as_int))
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
static void logger_bbob_open_file(FILE **file, const std::string &file_path) {
  if (*file == nullptr) {
    *file = fopen(file_path.c_str(), "a");
    if (*file == nullptr) {
      coco_error("logger_bbob_open_file(): Error opening file: %s\nError: %d", file_path.c_str(), errno);
    }
  }
}

/**
 * @brief Creates the data file (if it didn't exist before) and opens it
 */
static void logger_bbob_open_data_file(FILE **data_file, const std::string &path, const std::string &file_name,
                                       const char *file_extension) {
  logger_bbob_open_file(data_file, coco_join_path(path, file_name + file_extension));
}

/**
 * @brief Creates the info file (if it didn't exist before) and opens it
 */
void logger_bbob_t::open_info_file(const std::string &folder, const std::string &function_string,
                                   const std::string &data_file_name, int start_new_line) {
  int add_empty_line = 0;
  FILE *tmp_file;

  coco_debug("Started logger_bbob_open_info_file()");

  assert(observer != nullptr);

  const std::string file_path = coco_join_path(folder, observer->prefix + "_f" + function_string + ".info");

  if (info_file == nullptr) {
    add_empty_line = 0;
    /* If the file already exists, an empty line is needed */
    tmp_file = fopen(file_path.c_str(), "r");
    if (tmp_file != nullptr) {
      add_empty_line = 1;
      fclose(tmp_file);
    }
    logger_bbob_open_file(&info_file, file_path);
    if (start_new_line) {
      if (add_empty_line)
        fprintf(info_file, "\n");

      fprintf(info_file,
              "suite = '%s', funcId = %lu, DIM = %lu, Precision = %.3e, "
              "algId = '%s', coco_version = '%s', logger = '%s', "
              "data_format = '%s', settings = '%s'\n",
              suite_name.c_str(), (unsigned long)suite_dep_function, (unsigned long)number_of_variables, pow(10, -8),
              observer->algorithm_name.c_str(), coco_version, observer->observer_name.c_str(), logger_bbob_data_format,
              observer->settings.c_str());
      fprintf(info_file, "%% %s\n", observer->algorithm_info.c_str());
      /* data_file_name does not have the extension */
      fprintf(info_file, "%s.dat", data_file_name.c_str());
    }
  }
  coco_debug("Ended   logger_bbob_open_info_file()");
}

/**
 * @brief Generates the files and folders needed by the logger if they don't already exist
 */
void logger_bbob_t::initialize() {

  std::string relative_path; /* Relative path to the .dat file from where the .info file is */
  const char *str_pointer;
  int start_new_line;
  observer_bbob_t *observer_data;

  coco_debug("Started logger_bbob_initialize()");

  if (is_opt_known == 0)
    str_pointer = "Fref";
  else
    str_pointer = "Fopt";

  assert(observer != nullptr);
  observer_data = observer;

  const std::string function_string = coco_string_format("%lu", (unsigned long)suite_dep_function);
  const std::string dimension_string = coco_string_format("%lu", (unsigned long)number_of_variables);
  const std::string &result_folder = observer->result_folder;

  /* Prepare paths and names */
  coco_create_directory(coco_join_path(result_folder, "data_f" + function_string));
  relative_path = "data_f" + function_string + "/bbobexp_f" + function_string + "_DIM" + dimension_string;

  /* Find out whether a new line in the info file (and a new .?dat file) is needed */
  start_new_line = logger_bbob_start_new_line(observer, number_of_variables, suite_dep_function);
  if (start_new_line) {
    relative_path = coco_create_unique_filename(result_folder, relative_path, ".dat");
    observer_data->last_dat_file = relative_path;
    observer_data->last_function = suite_dep_function;
    observer_data->last_dimension = number_of_variables;
  } else {
    relative_path = observer_data->last_dat_file;
  }

  /* info file */
  open_info_file(result_folder, function_string, relative_path, start_new_line);
  fprintf(info_file, ", %lu", (unsigned long)suite_dep_instance);

  /* data files */
  logger_bbob_open_data_file(&(dat_file), result_folder, relative_path, ".dat");
  fprintf(dat_file, logger_bbob_header, str_pointer, optimal_value);
  logger_bbob_open_data_file(&(tdat_file), result_folder, relative_path, ".tdat");
  fprintf(tdat_file, logger_bbob_header, str_pointer, optimal_value);
  logger_bbob_open_data_file(&(rdat_file), result_folder, relative_path, ".rdat");
  fprintf(rdat_file, logger_bbob_header, str_pointer, optimal_value);
  logger_bbob_open_data_file(&(mdat_file), result_folder, relative_path, ".mdat");
  fprintf(mdat_file, logger_bbob_header, str_pointer, optimal_value);

  is_initialized = 1;

  coco_debug("Ended   logger_bbob_initialize()");
}

/**
 * @brief Evaluates the function, increases the number of evaluations and outputs information according to
 * observer options.
 */
void logger_bbob_t::evaluate_function(const double *x, double *y) {
  size_t i;
  double y_logged, max_value = 0, sum_constraints;
  const int feasible = number_of_constraints <= 0 || inner_problem->is_feasible(x, nullptr);

  if (is_tainted || inner_problem->is_tainted) {
    coco_error("Though shalt not cheat! Trying to evaluate a tainted observed problem. Do _not_ call "
               "coco_problem_best_parameter() on a problem you want to benchmark.");
  }

  coco_debug("Started logger_bbob_evaluate()");

  if (!is_initialized) {
    initialize();
  }
  if ((coco_log_level >= COCO_DEBUG) && num_func_evaluations == 0) {
    coco_debug("%4lu: ", (unsigned long)inner_problem->suite_dep_index);
    coco_debug("on problem %s ... ", coco_problem_get_id(inner_problem.get()));
  }

  /* Fulfill contract of a COCO evaluate function */
  coco_evaluate_function(inner_problem.get(), x, y);
  num_func_evaluations++;
  current_value = y[0];
  y_logged = y[0];
  if (inner_problem->is_noisy == 1) {
    current_value = inner_problem->last_noise_free_values[0];
    y_logged = inner_problem->last_noise_free_values[0];
  }
  if (coco_is_nan(y_logged))
    y_logged = NAN_FOR_LOGGING;
  else if (coco_is_inf(y_logged))
    y_logged = INFINITY_FOR_LOGGING;

  /* Do sanity check */
  if ((is_opt_known) && (feasible)) {
    /* Infeasible solutions can have much better y0 values */
    assert(y_logged + 1e-13 >= optimal_value);
  }

  /* Evaluate the constraints */
  std::vector<double> constraints(number_of_constraints);
  if (number_of_constraints > 0)
    inner_problem->evaluate_constraint(x, constraints.data(), 0);
  num_cons_evaluations = evaluations_constraints;

  /* Compute the sum of positive constraint values */
  sum_constraints = 0;
  for (i = 0; i < number_of_constraints; ++i) {
    if (constraints[i] > 0)
      sum_constraints += constraints[i];
  }
  sum_constraints *= LOGGER_BBOB_WEIGHT_CONSTRAINTS; /* do this before the checks */
  if (coco_is_nan(sum_constraints))
    sum_constraints = NAN_FOR_LOGGING;
  else if (coco_is_inf(sum_constraints))
    sum_constraints = INFINITY_FOR_LOGGING;

  if (is_opt_known)
    max_value = coco_double_max(y_logged, optimal_value);
  else {
    /* Problems with unknown optimal values */
    if (sum_constraints > 0) {
      max_value = LOGGER_BBOB_INFEASIBLE_PENALTY; /* Infeasible solution */
    } else
      max_value = y_logged; /* Feasible solution */
  }

  /* Update logger state
   * At number_of_evaluations == 1 the best_found_value is not initialized,
   * also compare to max_value to not potentially be thrown off by weird values in y[0]
   */
  if (num_func_evaluations == 1 || (max_value + sum_constraints < best_found_value)) {
    best_found_value = max_value + sum_constraints;
    for (i = 0; i < number_of_variables; i++)
      best_found_solution[i] = x[i]; /* May well be infeasible */

    /* Add a line in the .dat file is the performance target is reached by a feasible solution
     * and always at the first evaluation */
    if (num_func_evaluations == 1 || targets.trigger(best_found_value - optimal_value)) {
      output(dat_file, x, y_logged, best_found_value, constraints.data());
      last_logged_eval_dat = num_func_evaluations;
    }
  }

  /* Add a line in the .tdat file for specific number of evaluations*/
  if (evaluation_triggers.trigger(num_func_evaluations + num_cons_evaluations)) {
    output(tdat_file, x, y_logged, best_found_value, constraints.data());
    last_logged_eval_tdat = num_func_evaluations;
  }

  /* Add a line in the .rdat file if the algorithm was restarted */
  if (algorithm_restarted) {
    output(rdat_file, x, y_logged, best_found_value, constraints.data());
    algorithm_restarted = 0;
  }

  coco_debug("Ended   logger_bbob_evaluate()");
}

/**
 * @brief Evaluates the function and outputs information according to observer options to the file with
 * recommendations. The evaluation result is not returned and the evaluation counter is not increased.
 */
void logger_bbob_t::recommend_solution(const double *x) {
  double y_logged;
  const int feasible = number_of_constraints <= 0 || inner_problem->is_feasible(x, nullptr);

  if (!is_initialized) {
    initialize();
  }

  /* Fulfill contract of a COCO evaluate function, but do not increase the evaluation counters */
  std::vector<double> y(number_of_objectives);
  coco_evaluate_function(inner_problem.get(), x, y.data());

  y_logged = y[0];
  if (coco_is_nan(y_logged))
    y_logged = NAN_FOR_LOGGING;
  else if (coco_is_inf(y_logged))
    y_logged = INFINITY_FOR_LOGGING;

  /* Do sanity check */
  if ((is_opt_known) && (feasible)) {
    /* Infeasible solutions can have much better y0 values */
    assert(y_logged + 1e-13 >= optimal_value);
  }

  /* Evaluate the constraints */
  std::vector<double> constraints(number_of_constraints);
  if (number_of_constraints > 0)
    inner_problem->evaluate_constraint(x, constraints.data(), 0);

  /* Compute the sum of positive constraint values */
  double sum_constraints;
  sum_constraints = 0.0;
  for (size_t i = 0; i < number_of_constraints; ++i) {
    if (constraints[i] > 0)
      sum_constraints += constraints[i];
  }
  sum_constraints *= LOGGER_BBOB_WEIGHT_CONSTRAINTS; /* do this before the checks */
  if (coco_is_nan(sum_constraints))
    sum_constraints = NAN_FOR_LOGGING;
  else if (coco_is_inf(sum_constraints))
    sum_constraints = INFINITY_FOR_LOGGING;

  double max_value;
  if (is_opt_known)
    max_value = coco_double_max(y_logged, optimal_value);
  else {
    /* Problems with unknown optimal values */
    if (sum_constraints > 0) {
      max_value = LOGGER_BBOB_INFEASIBLE_PENALTY; /* Infeasible solution */
    } else
      max_value = y_logged; /* Feasible solution */
  }

  /* Update recommender state
   * The best_recommended_value is initialized to INFINITY_FOR_LOGGING.
   */
  if (max_value + sum_constraints < best_recommended_value) {
    for (size_t i = 0; i < number_of_variables; i++)
      best_recommended_solution[i] = x[i]; /* May well be infeasible */

    best_recommended_value = max_value + sum_constraints;
    /* Add a line in the .mdat file is the performance target is reached by a feasible solution */
    if (mdat_targets.trigger(best_recommended_value - optimal_value)) {
      output(mdat_file, x, y_logged, best_recommended_value, constraints.data());
      last_logged_eval_mdat = num_func_evaluations;
    }
  }
}

/**
 * @brief Takes care of some final output and closes the files.
 */
logger_bbob_t::~logger_bbob_t() {
  coco_debug("Started logger_bbob_t::~logger_bbob_t()");

  if ((coco_log_level >= COCO_DEBUG) && num_func_evaluations > 0) {
    coco_debug("best f=%e after %lu fevals (done observing)\n", best_found_value, (unsigned long)num_func_evaluations);
  }
  if (info_file != nullptr) {
    fprintf(info_file, ":%lu|%.1e", (unsigned long)num_func_evaluations, best_found_value - optimal_value);
    fclose(info_file);
    info_file = nullptr;
  }

  if (dat_file != nullptr) {
    if (num_func_evaluations != last_logged_eval_dat)
      output(dat_file, best_found_solution.data(), best_found_value, best_found_value, nullptr);
    fclose(dat_file);
    dat_file = nullptr;
  }

  if (tdat_file != nullptr) {
    if (num_func_evaluations != last_logged_eval_tdat)
      output(tdat_file, best_found_solution.data(), best_found_value, best_found_value, nullptr);
    fclose(tdat_file);
    tdat_file = nullptr;
  }

  if (rdat_file != nullptr) {
    fclose(rdat_file);
    rdat_file = nullptr;
  }

  if (mdat_file != nullptr) {
    if ((num_func_evaluations != last_logged_eval_mdat) && (last_logged_eval_mdat != 0))
      output(mdat_file, best_recommended_solution.data(), best_recommended_value, best_recommended_value, nullptr);
    fclose(mdat_file);
    mdat_file = nullptr;
  }

  /* The observer might not exist at this point */
  if ((observer != nullptr) && (observer->is_active == 1)) {
    /* If the observer still exists (if it does not, observed_problem does not matter any longer) */
    observer->observed_problem = nullptr;
  }

  coco_debug("Ended   logger_bbob_t::~logger_bbob_t()");
}

/**
 * @brief Saves the information that the algorithm has restarted
 */
void logger_bbob_t::signal_restart() {
  if (num_func_evaluations > 0)
    algorithm_restarted = 1;
}

/**
 * @brief Returns the bbob logger behind the given problem (which must be one).
 */
logger_bbob_t *logger_bbob_t::from(const coco_problem_t *problem) {
  assert(dynamic_cast<const logger_bbob_t *>(problem) != nullptr);
  return static_cast<logger_bbob_t *>(const_cast<coco_problem_t *>(problem));
}

void observer_bbob_t::signal_restart(coco_problem_t *problem) { logger_bbob_t::from(problem)->signal_restart(); }

/**
 * @brief Makes the observed logger (if any) forget this observer.
 */
observer_bbob_t::~observer_bbob_t() {
  coco_debug("Started observer_bbob_t::~observer_bbob_t()");
  /* Make sure that the observed problem's pointer to the observer points to NULL */
  if (observed_problem != nullptr) {
    observed_problem->observer = nullptr;
    observed_problem = nullptr;
  }
  coco_debug("Ended   observer_bbob_t::~observer_bbob_t()");
}

logger_bbob_t::logger_bbob_t(observer_bbob_t &observer, coco_problem_ptr_t inner)
    : coco_problem_transformed_t(std::move(inner), observer.observer_name) {
  coco_debug("Started logger_bbob()");

  assert(inner_problem.get());
  assert(inner_problem->suite);
  const coco_suite_t *problem_suite = inner_problem->suite;
  if (inner_problem->number_of_objectives != 1) {
    coco_warning("logger_bbob(): The bbob logger shouldn't be used to log a problem with %lu objectives",
                 (unsigned long)inner_problem->number_of_objectives);
  }
  if (observer.observed_problem != nullptr) {
    coco_error("logger_bbob(): The observed problem must be closed before a new problem can be observed");
  }
  this->observer = &observer;
  suite_name = inner_problem->suite->suite_name;
  is_initialized = 0;
  algorithm_restarted = 0;

  info_file = nullptr;
  dat_file = nullptr;
  tdat_file = nullptr;
  rdat_file = nullptr;
  mdat_file = nullptr;

  num_func_evaluations = 0;
  num_cons_evaluations = 0;
  last_logged_eval_mdat = 0;
  /* last_logged_evaluation = 0; */

  best_found_solution.resize(inner_problem->number_of_variables);
  best_found_value = DBL_MAX;
  optimal_value = inner_problem->best_value[0];

  current_value = DBL_MAX;

  best_recommended_solution.resize(inner_problem->number_of_variables);
  best_recommended_value = INFINITY_FOR_LOGGING;

  log_discrete_as_int = observer.log_discrete_as_int;

  /* Initialize triggers performance and evaluation triggers depending on the observer option */
  targets = coco_observer_targets_t(problem_suite->known_optima, observer.lin_target_precision,
                                    observer.number_target_triggers, observer.log_target_precision);
  evaluation_triggers = coco_observer_evaluations_t(
      observer.number_evaluation_triggers, observer.base_evaluation_triggers, inner_problem->number_of_variables);
  mdat_targets = coco_observer_targets_t(problem_suite->known_optima, observer.lin_target_precision,
                                         observer.number_target_triggers, observer.log_target_precision);

  coco_debug("Ended   logger_bbob()");

  observer.observed_problem = this;
}

coco_problem_ptr_t observer_bbob_t::add_logger(coco_problem_ptr_t problem) {
  return std::make_unique<logger_bbob_t>(*this, std::move(problem));
}
