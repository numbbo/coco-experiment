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

namespace coco {

/**
 * @brief Factor used in logged indicator (f-f*)^+ + sum_i g_i^+ in front of the sum
 */
static constexpr double LOGGER_BBOB_WEIGHT_CONSTRAINTS = 1e0;

/**
 * Penalty used when logging infeasible solutions of problems with unknown optimal value
 */
static constexpr double LOGGER_BBOB_INFEASIBLE_PENALTY = 1e10;

/**
 * Record variable values only when #variables is below this threshold, was: 22
 */
static constexpr size_t LOGGER_BBOB_NUMBER_OF_VARIABLES_THRESHOLD = 7;

/**
 * @brief The bbob logger, which wraps the observed problem.
 */
class LoggerBbob final : public Logger {
public:
  LoggerBbob(ObserverBbob &observer, ProblemPtr inner);

  ~LoggerBbob() override;

  void EvaluateFunction(std::span<const double> x, std::span<double> y) override;
  void RecommendSolution(std::span<const double> x) override;

  /** @brief Saves the information that the algorithm has restarted. */
  void SignalRestart() override;

  double best_found_value() const { return best_found_value_; }
  double optimal_value() const { return optimal_value_; }

  /** @brief Returns the triggers based on target values (which the tests use directly). */
  ObserverTargets &targets() { return targets_; }

  /** @brief Returns the observer (nullptr if it has been deleted). */
  ObserverBbob *observer() const { return static_cast<ObserverBbob *>(Logger::observer()); }

private:
  /** @brief Generates the files and folders needed by the logger if they don't already exist. */
  void Initialize();

  /** @brief Creates the info file (if it didn't exist before) and opens it. */
  void OpenInfoFile(const std::string &folder, const std::string &function_string, const std::string &data_file_name,
                    int start_new_line);

  /**
   * @brief Decides whether a new line in the info file (and therefore a new set of .?dat files) should be
   * started, which is the case when the current dimension or function differ from the last logged ones.
   */
  bool StartNewLine(size_t current_dimension, size_t current_function) const {
    assert(observer() != nullptr);
    return (current_function != observer()->last_function_) || (current_dimension != observer()->last_dimension_);
  }

  /** @brief Opens the file in append mode unless it is already open. */
  static void OpenFile(File &file, const std::string &file_path);

  /** @brief Opens (and creates if needed) the data file path/file_name + file_extension unless it is already open. */
  static void OpenDataFile(File &data_file, const std::string &path, const std::string &file_name,
                           const char *file_extension) {
    OpenFile(data_file, JoinPath(path, file_name + file_extension));
  }

  /** @brief Outputs a formatted line to a data file. */
  void Output(FILE *data_file, std::span<const double> x, double logged_value, double logged_best_value,
              std::span<const double> constraints) const;

  /** @brief The decision values are not transformed. */
  bool Inverse(std::span<const double> inner_x, std::span<double> x) const override {
    return IdentityInverse(inner_x, x);
  }

  std::string suite_name_;        /**< @brief The suite name */
  bool is_initialized_{};         /**< @brief Whether the logger was already initialized */
  bool algorithm_restarted_{};    /**< @brief Whether the algorithm has restarted (output information to .rdat file). */
  File info_file_;                /**< @brief Index file */
  File dat_file_;                 /**< @brief File with function value aligned data */
  File tdat_file_;                /**< @brief File with number of evaluations aligned data */
  File rdat_file_;                /**< @brief File with restart information */
  File mdat_file_;                /**< @brief File with evaluated recommendations */
  size_t num_func_evaluations_{}; /**< @brief The number of function evaluations performed so far. */
  size_t num_cons_evaluations_{}; /**< @brief The number of evaluations of constraints performed so far. */
  size_t last_logged_eval_dat_{}; /**< @brief The last evaluation logged in the dat file (needed for finalization) */
  size_t last_logged_eval_tdat_{}; /**< @brief The last evaluation logged in the tdat file (needed for finalization) */
  size_t last_logged_eval_mdat_{}; /**< @brief The last evaluation logged in the mdat file (needed for finalization) */
  std::vector<double> best_found_solution_; /**< @brief The best found solution to this problem. */
  double best_found_value_{};               /**< @brief The best found value for this problem. */
  double current_value_{};                  /**< @brief The current value for this problem. */
  double
      optimal_value_{}; /**< @brief The optimal value for this problem (if it is known, otherwise a reference value). */
  std::vector<double> best_recommended_solution_; /**< @brief The best recommended solution for this problem. */
  double best_recommended_value_{};               /**< @brief The best recommended value for this problem. */
  bool log_discrete_as_int_{};      /**< @brief Whether to output discrete variables in integer or double format. */
  ObserverTargets targets_;         /**< @brief Triggers based on target values. */
  ObserverEvaluations evaluation_triggers_; /**< @brief Triggers based on the number of evaluations. */
  ObserverTargets mdat_targets_;            /**< @brief Triggers based on target values for recommendations. */
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
static const char *LOGGER_BBOB_DATA_FORMAT = "bbob-new2";

/**
 * @brief The header used by the logger in .?dat files
 *
 * See also logger_bbob_data_format
 */
static const char *LOGGER_BBOB_HEADER = "%% "
                                        "f evaluations | "
                                        "g evaluations | "
                                        "best noise-free fitness - %s (%13.12e) + sum g_i+ | "
                                        "measured fitness | "
                                        "best measured fitness or single-digit g-values | "
                                        "x1 | "
                                        "x2...\n";

void LoggerBbob::OpenFile(File &file, const std::string &file_path) {
  if (file == nullptr) {
    if (!file.Open(file_path, "a")) {
      coco_error("LoggerBbob::OpenFile(): Error opening file: %s\nError: %d", file_path.c_str(), errno);
    }
  }
}

/**
 * @brief Outputs a formatted line to a data file
 */
void LoggerBbob::Output(FILE *data_file, std::span<const double> x, double logged_value, double logged_best_value,
                        std::span<const double> constraints) const {
  /* This function contains many hard-coded values (10.9, 22, 5.4) that could be read through
   * observer options */

  fprintf(data_file, "%zu %zu %+10.9e %+10.9e ", num_func_evaluations_, num_cons_evaluations_,
          logged_best_value - optimal_value_, logged_value);

  if ((number_of_constraints_ > 0) && !constraints.empty()) {
    for (size_t i = 0; i < number_of_constraints_; ++i) {
      /* print 01234567890123..., may happen in the last line of .tdat */
      fprintf(data_file, "%d",
              !constraints.empty() ? SingleDigitConstraintValue(constraints[i]) : static_cast<int>(i % 10));
    }
  } else {
    fprintf(data_file, "%+10.9e", logged_best_value);
  }

  if (number_of_variables_ < LOGGER_BBOB_NUMBER_OF_VARIABLES_THRESHOLD) {
    for (size_t i = 0; i < number_of_variables_; i++) {
      if ((i < number_of_integer_variables_) && (log_discrete_as_int_))
        fprintf(data_file, " %d", DoubleToInt(x[i]));
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
 * @brief Creates the info file (if it didn't exist before) and opens it
 */
void LoggerBbob::OpenInfoFile(const std::string &folder, const std::string &function_string,
                              const std::string &data_file_name, int start_new_line) {
  int add_empty_line = 0;

  coco_debug("Started logger_bbob_t::OpenInfoFile()");

  assert(observer() != nullptr);

  const std::string file_path = JoinPath(folder, observer()->prefix_ + "_f" + function_string + ".info");

  if (info_file_ == nullptr) {
    add_empty_line = 0;
    /* If the file already exists, an empty line is needed */
    File tmp_file;
    if (tmp_file.Open(file_path, "r")) {
      add_empty_line = 1;
    }
    OpenFile(info_file_, file_path);
    if (start_new_line) {
      if (add_empty_line)
        fprintf(info_file_, "\n");

      fprintf(info_file_,
              "suite = '%s', funcId = %zu, DIM = %zu, Precision = %.3e, "
              "algId = '%s', coco_version = '%s', logger = '%s', "
              "data_format = '%s', settings = '%s'\n",
              suite_name_.c_str(), suite_dep_function(), number_of_variables_, pow(10, -8),
              observer()->algorithm_name.c_str(), coco_version, observer()->observer_name.c_str(),
              LOGGER_BBOB_DATA_FORMAT, observer()->settings.c_str());
      fprintf(info_file_, "%% %s\n", observer()->algorithm_info.c_str());
      /* data_file_name does not have the extension */
      fprintf(info_file_, "%s.dat", data_file_name.c_str());
    }
  }
  coco_debug("Ended   logger_bbob_t::OpenInfoFile()");
}

/**
 * @brief Generates the files and folders needed by the logger if they don't already exist
 */
void LoggerBbob::Initialize() {

  std::string relative_path; /* Relative path to the .dat file from where the .info file is */
  const char *str_pointer;
  ObserverBbob *observer_data;

  coco_debug("Started logger_bbob_t::initialize()");

  if (opt_known_ == 0)
    str_pointer = "Fref";
  else
    str_pointer = "Fopt";

  assert(observer() != nullptr);
  observer_data = observer();

  const std::string function_string = StringFormat("%zu", suite_dep_function());
  const std::string dimension_string = StringFormat("%zu", number_of_variables_);
  const std::string &result_folder = observer()->result_folder;

  /* Prepare paths and names */
  CreateDirectory(JoinPath(result_folder, "data_f" + function_string));
  relative_path = "data_f" + function_string + "/bbobexp_f" + function_string + "_DIM" + dimension_string;

  /* Find out whether a new line in the info file (and a new .?dat file) is needed */
  const bool new_line = StartNewLine(number_of_variables_, suite_dep_function());
  if (new_line) {
    relative_path = CreateUniqueFilename(result_folder, relative_path, ".dat");
    observer_data->last_dat_file_ = relative_path;
    observer_data->last_function_ = suite_dep_function();
    observer_data->last_dimension_ = number_of_variables_;
  } else {
    relative_path = observer_data->last_dat_file_;
  }

  /* info file */
  OpenInfoFile(result_folder, function_string, relative_path, new_line);
  fprintf(info_file_, ", %zu", suite_dep_instance());

  /* data files */
  OpenDataFile(dat_file_, result_folder, relative_path, ".dat");
  fprintf(dat_file_, LOGGER_BBOB_HEADER, str_pointer, optimal_value_);
  OpenDataFile(tdat_file_, result_folder, relative_path, ".tdat");
  fprintf(tdat_file_, LOGGER_BBOB_HEADER, str_pointer, optimal_value_);
  OpenDataFile(rdat_file_, result_folder, relative_path, ".rdat");
  fprintf(rdat_file_, LOGGER_BBOB_HEADER, str_pointer, optimal_value_);
  OpenDataFile(mdat_file_, result_folder, relative_path, ".mdat");
  fprintf(mdat_file_, LOGGER_BBOB_HEADER, str_pointer, optimal_value_);

  is_initialized_ = 1;

  coco_debug("Ended   logger_bbob_t::initialize()");
}

/**
 * @brief Evaluates the function, increases the number of evaluations and outputs information according to
 * observer options.
 */
void LoggerBbob::EvaluateFunction(std::span<const double> x, std::span<double> y) {
  double y_logged, max_value = 0, sum_constraints;
  const bool feasible = number_of_constraints_ <= 0 || inner_problem_->IsFeasible(x, {});

  if (is_tainted() || inner_problem_->is_tainted()) {
    coco_error("Though shalt not cheat! Trying to evaluate a tainted observed problem. Do _not_ call "
               "coco_problem_best_parameter() on a problem you want to benchmark.");
  }

  coco_debug("Started logger_bbob_t::EvaluateFunction()");

  if (!is_initialized_) {
    Initialize();
  }
  if ((coco_log_level >= COCO_DEBUG) && num_func_evaluations_ == 0) {
    coco_debug("%4zu: ", inner_problem_->suite_dep_index());
    coco_debug("on problem %s ... ", coco_problem_get_id(inner_problem_.get()));
  }

  /* Fulfill contract of a COCO evaluate function */
  inner_problem_->EvaluateFunction(x, y);
  num_func_evaluations_++;
  current_value_ = y[0];
  y_logged = y[0];
  if (inner_problem_->is_noisy() == 1) {
    current_value_ = inner_problem_->last_noise_free_values()[0];
    y_logged = inner_problem_->last_noise_free_values()[0];
  }
  if (IsNan(y_logged))
    y_logged = NAN_FOR_LOGGING;
  else if (IsInf(y_logged))
    y_logged = INFINITY_FOR_LOGGING;

  /* Do sanity check */
  if ((opt_known_) && (feasible)) {
    /* Infeasible solutions can have much better y0 values */
    assert(y_logged + 1e-13 >= optimal_value_);
  }

  /* Evaluate the constraints */
  std::vector<double> constraints(number_of_constraints_);
  if (number_of_constraints_ > 0)
    inner_problem_->EvaluateConstraint(x, constraints, 0);
  num_cons_evaluations_ = evaluations_constraints();

  /* Compute the sum of positive constraint values */
  sum_constraints = 0;
  for (size_t i = 0; i < number_of_constraints_; ++i) {
    if (constraints[i] > 0)
      sum_constraints += constraints[i];
  }
  sum_constraints *= LOGGER_BBOB_WEIGHT_CONSTRAINTS; /* do this before the checks */
  if (IsNan(sum_constraints))
    sum_constraints = NAN_FOR_LOGGING;
  else if (IsInf(sum_constraints))
    sum_constraints = INFINITY_FOR_LOGGING;

  if (opt_known_)
    max_value = DoubleMax(y_logged, optimal_value_);
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
  if (num_func_evaluations_ == 1 || (max_value + sum_constraints < best_found_value_)) {
    best_found_value_ = max_value + sum_constraints;
    for (size_t i = 0; i < number_of_variables_; i++)
      best_found_solution_[i] = x[i]; /* May well be infeasible */

    /* Add a line in the .dat file is the performance target is reached by a feasible solution
     * and always at the first evaluation */
    if (num_func_evaluations_ == 1 || targets_.Trigger(best_found_value_ - optimal_value_)) {
      Output(dat_file_, x, y_logged, best_found_value_, constraints);
      last_logged_eval_dat_ = num_func_evaluations_;
    }
  }

  /* Add a line in the .tdat file for specific number of evaluations*/
  if (evaluation_triggers_.Trigger(num_func_evaluations_ + num_cons_evaluations_)) {
    Output(tdat_file_, x, y_logged, best_found_value_, constraints);
    last_logged_eval_tdat_ = num_func_evaluations_;
  }

  /* Add a line in the .rdat file if the algorithm was restarted */
  if (algorithm_restarted_) {
    Output(rdat_file_, x, y_logged, best_found_value_, constraints);
    algorithm_restarted_ = 0;
  }

  coco_debug("Ended   logger_bbob_t::EvaluateFunction()");
}

/**
 * @brief Evaluates the function and outputs information according to observer options to the file with
 * recommendations. The evaluation result is not returned and the evaluation counter is not increased.
 */
void LoggerBbob::RecommendSolution(std::span<const double> x) {
  const bool feasible = number_of_constraints_ <= 0 || inner_problem_->IsFeasible(x, {});

  if (!is_initialized_) {
    Initialize();
  }

  /* Fulfill contract of a COCO evaluate function, but do not increase the evaluation counters */
  std::vector<double> y(number_of_objectives_);
  coco_evaluate_function(inner_problem_.get(), x.data(), y.data()); /* checks x, which recommendations do not */

  double y_logged = y[0];
  if (IsNan(y_logged))
    y_logged = NAN_FOR_LOGGING;
  else if (IsInf(y_logged))
    y_logged = INFINITY_FOR_LOGGING;

  /* Do sanity check */
  if ((opt_known_) && (feasible)) {
    /* Infeasible solutions can have much better y0 values */
    assert(y_logged + 1e-13 >= optimal_value_);
  }

  /* Evaluate the constraints */
  std::vector<double> constraints(number_of_constraints_);
  if (number_of_constraints_ > 0)
    inner_problem_->EvaluateConstraint(x, constraints, 0);

  /* Compute the sum of positive constraint values */
  double sum_constraints = 0.0;
  for (size_t i = 0; i < number_of_constraints_; ++i) {
    if (constraints[i] > 0)
      sum_constraints += constraints[i];
  }
  sum_constraints *= LOGGER_BBOB_WEIGHT_CONSTRAINTS; /* do this before the checks */
  if (IsNan(sum_constraints))
    sum_constraints = NAN_FOR_LOGGING;
  else if (IsInf(sum_constraints))
    sum_constraints = INFINITY_FOR_LOGGING;

  double max_value;
  if (opt_known_)
    max_value = DoubleMax(y_logged, optimal_value_);
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
  if (max_value + sum_constraints < best_recommended_value_) {
    for (size_t i = 0; i < number_of_variables_; i++)
      best_recommended_solution_[i] = x[i]; /* May well be infeasible */

    best_recommended_value_ = max_value + sum_constraints;
    /* Add a line in the .mdat file is the performance target is reached by a feasible solution */
    if (mdat_targets_.Trigger(best_recommended_value_ - optimal_value_)) {
      Output(mdat_file_, x, y_logged, best_recommended_value_, constraints);
      last_logged_eval_mdat_ = num_func_evaluations_;
    }
  }
}

/**
 * @brief Takes care of some final output and closes the files.
 */
LoggerBbob::~LoggerBbob() {
  coco_debug("Started logger_bbob_t::~LoggerBbob()");

  if ((coco_log_level >= COCO_DEBUG) && num_func_evaluations_ > 0) {
    coco_debug("best f=%e after %zu fevals (done observing)\n", best_found_value_, num_func_evaluations_);
  }
  if (info_file_ != nullptr) {
    fprintf(info_file_, ":%zu|%.1e", num_func_evaluations_, best_found_value_ - optimal_value_);
    info_file_.Close();
  }

  if (dat_file_ != nullptr) {
    if (num_func_evaluations_ != last_logged_eval_dat_)
      Output(dat_file_, best_found_solution_, best_found_value_, best_found_value_, {});
    dat_file_.Close();
  }

  if (tdat_file_ != nullptr) {
    if (num_func_evaluations_ != last_logged_eval_tdat_)
      Output(tdat_file_, best_found_solution_, best_found_value_, best_found_value_, {});
    tdat_file_.Close();
  }

  if (mdat_file_ != nullptr) {
    if ((num_func_evaluations_ != last_logged_eval_mdat_) && (last_logged_eval_mdat_ != 0))
      Output(mdat_file_, best_recommended_solution_, best_recommended_value_, best_recommended_value_, {});
    mdat_file_.Close();
  }

  coco_debug("Ended   logger_bbob_t::~LoggerBbob()");
}

/**
 * @brief Saves the information that the algorithm has restarted
 */
void LoggerBbob::SignalRestart() {
  if (num_func_evaluations_ > 0)
    algorithm_restarted_ = 1;
}

LoggerBbob::LoggerBbob(ObserverBbob &observer, ProblemPtr inner) : Logger(observer, std::move(inner)) {
  coco_debug("Started LoggerBbob()");

  assert(inner_problem_.get());
  assert(inner_problem_->suite());
  const Suite *problem_suite = inner_problem_->suite();
  if (inner_problem_->number_of_objectives() != 1) {
    coco_warning("LoggerBbob(): The bbob logger shouldn't be used to log a problem with %zu objectives",
                 inner_problem_->number_of_objectives());
  }
  suite_name_ = inner_problem_->suite()->name();
  is_initialized_ = 0;
  algorithm_restarted_ = 0;

  num_func_evaluations_ = 0;
  num_cons_evaluations_ = 0;
  last_logged_eval_mdat_ = 0;

  best_found_solution_.resize(inner_problem_->number_of_variables());
  best_found_value_ = DBL_MAX;
  optimal_value_ = inner_problem_->best_value()[0];

  current_value_ = DBL_MAX;

  best_recommended_solution_.resize(inner_problem_->number_of_variables());
  best_recommended_value_ = INFINITY_FOR_LOGGING;

  log_discrete_as_int_ = observer.log_discrete_as_int;

  /* Initialize triggers performance and evaluation triggers depending on the observer option */
  targets_ = ObserverTargets(problem_suite->has_known_optima(), observer.lin_target_precision,
                             observer.number_target_triggers, observer.log_target_precision);
  evaluation_triggers_ = ObserverEvaluations(observer.number_evaluation_triggers, observer.base_evaluation_triggers,
                                             inner_problem_->number_of_variables());
  mdat_targets_ = ObserverTargets(problem_suite->has_known_optima(), observer.lin_target_precision,
                                  observer.number_target_triggers, observer.log_target_precision);

  coco_debug("Ended   LoggerBbob()");
}

ProblemPtr ObserverBbob::AddLogger(ProblemPtr problem) {
  if (number_of_loggers() > 0)
    coco_error("LoggerBbob(): The observed problem must be closed before a new problem can be observed");
  return std::make_unique<LoggerBbob>(*this, std::move(problem));
}

} // namespace coco
