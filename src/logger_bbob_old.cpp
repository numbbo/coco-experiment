/**
 * @file logger_bbob_old.cpp
 * @brief Old implementation of the bbob logger.
 *
 * Logs the performance of a single-objective optimizer on noisy or noiseless problems.
 * It produces four kinds of files:
 * - The "info" files ...
 * - The "dat" files ...
 * - The "tdat" files ...
 * - The "rdat" files ...
 */

/* TODO: Document this file in doxygen style! */

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
#include "observer_bbob_old.cpp"

namespace coco {

static const double FVALUE_LOGGED_FOR_INFINITE = 3e21; /* value used for logging try */
static const double FVALUE_LOGGED_FOR_NAN = 2e21;
static const double WEIGHT_CONSTRAINTS =
    1e0; /* factor used in logged indicator (f-f*)^+ + sum_i g_i^+ in front of the sum */

/* TODO: add possibility of adding a prefix to the index files (easy to do through observer options) */

/**
 * @brief The bbob-old logger, which wraps the observed problem.
 */
class LoggerBbobOld final : public Logger {
public:
  /*char *path;// relative path to the data folder. //Wassim: now fetched from the observer */

  LoggerBbobOld(ObserverBbobOld &observer, ProblemPtr inner);

  ~LoggerBbobOld() override;

  void EvaluateFunction(std::span<const double> x, std::span<double> y) override;

  /** @brief The decision values are not transformed. */
  bool Inverse(std::span<const double> inner_x, std::span<double> x) const override {
    return IdentityInverse(inner_x, x);
  }

private:
  /** @brief Generates the files and folders needed by the logger if they don't already exist. */
  void Initialize();

  /** @brief Opens the file in append mode, reports an error if that fails. */
  static void OpenFile(File &file, const std::string &file_path);

  /** @brief Opens (and creates if needed) the data file result_folder_/data_file_path + file_extension. */
  void OpenDataFile(File &file, const std::string &data_file_path, const char *file_extension) {
    if (file == nullptr)
      OpenFile(file, JoinPath(result_folder_, data_file_path + file_extension));
  }

  /** @brief Adds a formatted line to the data file. */
  void WriteData(FILE *file, double fvalue, double best_fvalue, std::span<const double> x,
                 std::span<const double> constraints) const;

  /** @brief Creates the index file fileName_prefix+problem_id+file_extension in folder_path. */
  void OpenIndexFile(const std::string &folder_path, const char *indexFile_prefix, const std::string &function_id,
                     const std::string &dataFile_path, const std::string &suite_name);

  std::shared_ptr<LoggerBbobOldInfoState> info_state_; /**< @brief Shared with the observer */
  std::string result_folder_;  /**< @brief Copied from the observer, which might be freed before the logger. */
  std::string algorithm_name_; /**< @brief Copied from the observer, which might be freed before the logger. */
  bool is_initialized_{};
  /*const char *alg_name; the alg name, for now, temporarily the same as the path. Wassim: Now in the observer */
  File index_file_; /* index file */
  File fdata_file_; /* function value aligned data file */
  File tdata_file_; /* number of function evaluations aligned data file */
  File rdata_file_; /* restart info data file */
  size_t number_of_evaluations_{};
  size_t number_of_evaluations_constraints_{};
  double best_fvalue_{};
  double last_fvalue_{};
  bool written_last_eval_{}; /* allows writing the data of the final fun eval in the .tdat file if not already written
                              by the t_trigger*/
  std::vector<double> best_solution_;
  bool log_discrete_as_int_{}; /**< @brief Whether to output discrete variables in int or double format. */
  double optimal_fvalue_{};
  ObserverLogTargets targets_;              /**< @brief Triggers based on target values. */
  ObserverEvaluations evaluation_triggers_; /**< @brief Triggers based on the number of evaluations. */
};

/* was (old):
 * function evaluation |
 * noise-free fitness - Fopt (7.948000000000e+01) |
 * best noise-free fitness - Fopt |
 * measured fitness |
 * best measured fitness |
 * x1 | x2...
 was (bbob-new):
    "f evaluations | "
    "g evaluations | "
    "best noise-free fitness - Fopt | "
    "noise-free fitness - Fopt (%13.12e) | "
    "measured fitness | "
    "best measured fitness | "
    "x1 | "
    "x2...\n";
 */
static const char *BBOB_OLD_FILE_HEADER_STR = "%% "
                                              "f evaluations | "
                                              "g evaluations | "
                                              "best noise-free fitness - Fopt (%13.12e) + sum g_i+ | "
                                              "measured fitness | "
                                              "best measured fitness or single-digit g-values | "
                                              "x1 | "
                                              "x2...\n";

static const char *LOGGER_NAME = "bbob_old";
static const char *DATA_FORMAT = "bbob-new2"; /* back to 5 columns, 5-th column writes single digit constraint values */

void LoggerBbobOld::OpenFile(File &file, const std::string &file_path) {
  if (!file.Open(file_path, "a+"))
    coco_error("LoggerBbobOld::OpenFile(): Error opening file %s: %s", file_path.c_str(), strerror(errno));
}

void LoggerBbobOld::WriteData(FILE *file, double fvalue, double best_fvalue, std::span<const double> x,
                              std::span<const double> constraints) const {
  /* for some reason, it's %.0f in the old code instead of the 10.9e
   * in the documentation
   */
  fprintf(file, "%zu %zu %+10.9e %+10.9e ", number_of_evaluations_, number_of_evaluations_constraints_,
          best_fvalue - optimal_fvalue_, fvalue);

  if (!constraints.empty())
    for (const double constraint : constraints)
      fprintf(file, "%d", SingleDigitConstraintValue(constraint));
  else
    fprintf(file, "%+10.9e", best_fvalue);

  if ((x.size() - number_of_integer_variables_) < 22) {
    for (size_t i = 0; i < x.size(); i++) {
      if ((i < number_of_integer_variables_) && (log_discrete_as_int_))
        fprintf(file, " %d", DoubleToInt(x[i]));
      else
        fprintf(file, " %+5.4e", x[i]);
    }
  }
  fprintf(file, "\n");

  /* Flush output so that impatient users can see progress.
   * Otherwise it can take a long time until the output appears.
   */
  fflush(file);
}

void LoggerBbobOld::OpenIndexFile(const std::string &folder_path, const char *indexFile_prefix,
                                  const std::string &function_id, const std::string &dataFile_path,
                                  const std::string &suite_name) {
  /* to add the instance number TODO: this should be done outside to avoid redoing this for the .*dat files */
  std::string used_dataFile_path = dataFile_path;
  int newLine = 0;         /* newLine is at 1 if we need a new line in the info file */
  File tmp_file;           /* to check whether the file already exists */
  if (info_state_->first_instance == 0) {
    info_state_->first_instance = suite_dep_instance();
  }
  const std::string function_id_char = StringFormat("%zu", suite_dep_function());
  info_state_->first_instance_string = StringFormat("%zu", info_state_->first_instance);
  File &target_file = index_file_;
  std::string file_path = JoinPath(folder_path, std::string(indexFile_prefix) + "_f" + function_id_char + "_i" +
                                                    info_state_->first_instance_string + ".info");
  if (target_file == nullptr) {
    tmp_file.Open(file_path, "r"); /* to check for existence */
    if ((tmp_file) && (info_state_->current_dimension == number_of_variables_) &&
        (info_state_->current_function == suite_dep_function())) {
      /* new instance of current funId and current dim */
      OpenFile(target_file, file_path);
      tmp_file.Close();
    } else { /* either file doesn't exist (new funId) or new Dim */
      /* check that the dim was not already present earlier in the file, if so, create a new info file */
      if (info_state_->current_dimension != number_of_variables_) {
        size_t i;
        for (i = 0; i < LoggerBbobOldInfoState::NUMBER_OF_DIMENSIONS && info_state_->dimensions_in_info_file[i] != 0 &&
                    info_state_->dimensions_in_info_file[i] != number_of_variables_;
             i++) {
          ; /* checks whether dimension already present in the current infoFile */
        }
        if (i < LoggerBbobOldInfoState::NUMBER_OF_DIMENSIONS && info_state_->dimensions_in_info_file[i] == 0) {
          /* new dimension seen for the first time */
          info_state_->dimensions_in_info_file[i] = number_of_variables_;
          newLine = 1;
        } else {
          if (i < LoggerBbobOldInfoState::NUMBER_OF_DIMENSIONS) { /* dimension already present, need to create a
                                                                           new file */
            newLine = 0;
            /* truncate the instance part ("_i" + instance + ".info") */
            file_path.resize(file_path.size() - info_state_->first_instance_string.size() - 7);
            info_state_->first_instance = suite_dep_instance();
            info_state_->first_instance_string = StringFormat("%zu", info_state_->first_instance);
            file_path += "_i" + info_state_->first_instance_string + ".info";
          } else { /* we have all dimensions */
            newLine = 1;
          }
          for (size_t j = 0; j < LoggerBbobOldInfoState::NUMBER_OF_DIMENSIONS;
               j++) { /* new info file, reinitialize list of dims */
            info_state_->dimensions_in_info_file[j] = 0;
          }
          info_state_->dimensions_in_info_file[i] = number_of_variables_;
        }
      } else {
        if (info_state_->current_function != suite_dep_function()) {
          /*new function in the same file */
          newLine = 1;
        }
      }
      OpenFile(target_file, file_path); /* in any case, we append */
      if (tmp_file) { /* File already exists, new dim so just a new line. Also, close the tmp_file */
        if (newLine) {
          fprintf(target_file, "\n");
        }
        tmp_file.Close();
      }
      fprintf(target_file,
              "suite = '%s', funcId = %d, DIM = %zu, Precision = %.3e, algId = '%s', coco_version = '%s', logger = "
              "'%s', data_format = '%s'\n",
              suite_name.c_str(), static_cast<int>(strtol(function_id.c_str(), nullptr, 10)), number_of_variables_,
              pow(10, -8), algorithm_name_.c_str(), coco_version, LOGGER_NAME, DATA_FORMAT);

      fprintf(target_file, "%%\n");
      used_dataFile_path += "_i" + info_state_->first_instance_string;
      fprintf(target_file, "%s.dat", used_dataFile_path.c_str()); /* dataFile_path does not have the extension */
      info_state_->current_dimension = number_of_variables_;
      info_state_->current_function = suite_dep_function();
    }
  }
}

void LoggerBbobOld::Initialize() {
  /*
   Creates/opens the data and index files
   */
  std::string dataFile_path; /* relative path to the .dat file from where the .info file is */
  char indexFile_prefix[10] =
      "bbobexp"; /* TODO (minor): make the prefix bbobexp a parameter that the user can modify */

  assert(inner_problem_ != nullptr);
  assert(!inner_problem_->id().empty());

  const std::string tmpc_funId = StringFormat("%zu", inner_problem_->suite_dep_function());
  const std::string tmpc_dim = StringFormat("%zu", inner_problem_->number_of_variables());

  /* prepare paths and names */
  CreateDirectory(JoinPath(result_folder_, "data_f" + tmpc_funId));
  dataFile_path = "data_f" + tmpc_funId + "/bbobexp_f" + tmpc_funId + "_DIM" + tmpc_dim;

  /* index/info file */
  assert(inner_problem_->suite());
  OpenIndexFile(result_folder_, indexFile_prefix, tmpc_funId, dataFile_path, inner_problem_->suite()->name());
  fprintf(index_file_, ", %zu", inner_problem_->suite_dep_instance());
  /* data files */
  /* TODO: definitely improvable but works for now */
  dataFile_path += "_i" + info_state_->first_instance_string;

  OpenDataFile(fdata_file_, dataFile_path, ".dat");
  fprintf(fdata_file_, BBOB_OLD_FILE_HEADER_STR, optimal_fvalue_);

  OpenDataFile(tdata_file_, dataFile_path, ".tdat");
  fprintf(tdata_file_, BBOB_OLD_FILE_HEADER_STR, optimal_fvalue_);

  OpenDataFile(rdata_file_, dataFile_path, ".rdat");
  fprintf(rdata_file_, BBOB_OLD_FILE_HEADER_STR, optimal_fvalue_);
  is_initialized_ = true;
}

/**
 * Layer added to the transformed-problem EvaluateFunction by the logger
 */
void LoggerBbobOld::EvaluateFunction(std::span<const double> x, std::span<double> y) {
  double y_logged, max_fvalue, sum_cons;
  const bool feasible = number_of_constraints_ <= 0 || inner_problem_->IsFeasible(x, {});

  if (!is_initialized_) {
    Initialize();
  }
  if ((coco_log_level >= COCO_DEBUG) && number_of_evaluations_ == 0) {
    coco_debug("%4zu: ", inner_problem_->suite_dep_index());
    coco_debug("on problem %s ... ", coco_problem_get_id(inner_problem_.get()));
  }

  inner_problem_->EvaluateFunction(x, y); /* fulfill contract as "being" a coco evaluate function */

  number_of_evaluations_constraints_ = coco_problem_get_evaluations_constraints(this);
  number_of_evaluations_++; /* could be != coco_problem_get_evaluations(problem) for non-anytime logging? */
  written_last_eval_ = 0;   /* flag whether the current evaluation was logged? */
  last_fvalue_ = y[0];      /* asma: should be: max(y[0], optimal_fvalue) */

  y_logged = y[0];
  if (IsNan(y_logged))
    y_logged = FVALUE_LOGGED_FOR_NAN;
  else if (IsInf(y_logged))
    y_logged = FVALUE_LOGGED_FOR_INFINITE;
  /* do sanity check */
  if (feasible) /* infeasible solutions can have much better y0 values */
    assert(y_logged + 1e-13 >= optimal_fvalue_);

  /* Evaluate the constraints */
  std::vector<double> cons(number_of_constraints_);
  if (number_of_constraints_ > 0)
    inner_problem_->EvaluateConstraint(x, cons, 0);

  /* Compute the sum of positive constraint values */
  sum_cons = 0;
  for (size_t i = 0; i < number_of_constraints_; ++i) {
    if (cons[i] > 0)
      sum_cons += cons[i];
  }
  sum_cons *= WEIGHT_CONSTRAINTS; /* do this before the checks */
  if (IsNan(sum_cons))
    sum_cons = FVALUE_LOGGED_FOR_NAN;
  else if (IsInf(sum_cons))
    sum_cons = FVALUE_LOGGED_FOR_INFINITE;

  max_fvalue = y_logged > optimal_fvalue_ ? y_logged : optimal_fvalue_;

  /* Update logger state.
   *   at number_of_evaluations == 1 the best_fvalue is not initialized,
   *   also compare to y_logged to not potentially be thrown off by weird values in y[0]
   */
  if (number_of_evaluations_ == 1 || (max_fvalue + sum_cons < best_fvalue_)) {
    best_fvalue_ = max_fvalue + sum_cons;
    for (size_t i = 0; i < number_of_variables_; i++)
      best_solution_[i] = x[i]; /* may well be infeasible */

    /* Add a line in the .dat file for each logging target reached
     * by a feasible solution and always at evaluation one
     */
    if (number_of_evaluations_ == 1 || targets_.Trigger(best_fvalue_ - optimal_fvalue_)) {
      WriteData(fdata_file_, y_logged, best_fvalue_, x, cons);
    }
  }

  /* Add a line in the .tdat file each time an fevals trigger is reached.*/
  if (evaluation_triggers_.Trigger(number_of_evaluations_ + number_of_evaluations_constraints_)) {
    WriteData(tdata_file_, y_logged, best_fvalue_, x, cons);
    written_last_eval_ = 1;
  }

} /* end logger_bbob_old_evaluate */

/**
 * @brief Writes the final output and closes the files.
 */
LoggerBbobOld::~LoggerBbobOld() {
  if ((coco_log_level >= COCO_DEBUG) && number_of_evaluations_ > 0) {
    coco_debug("best f=%e after %zu fevals (done observing)\n", best_fvalue_, number_of_evaluations_);
  }
  if (index_file_ != nullptr) {
    fprintf(index_file_, ":%zu|%.1e", number_of_evaluations_, best_fvalue_ - optimal_fvalue_);
    index_file_.Close();
  }
  if (tdata_file_ != nullptr) {
    /* TODO: make sure it handles restarts well. i.e., it writes
     * at the end of a single run, not all the runs on a given
     * instance. Maybe start with forcing it to generate a new
     * "instance" of problem for each restart in the beginning
     */
    if (!written_last_eval_) {
      WriteData(tdata_file_, best_fvalue_, best_fvalue_, best_solution_, {});
    }
    tdata_file_.Close();
  }
}

LoggerBbobOld::LoggerBbobOld(ObserverBbobOld &observer, ProblemPtr inner) : Logger(observer, std::move(inner)) {
  result_folder_ = observer.result_folder;
  info_state_ = observer.info_state_;
  algorithm_name_ = observer.algorithm_name;

  if (inner_problem_->number_of_objectives() != 1) {
    coco_warning("LoggerBbobOld(): The bbob_old logger shouldn't be used to log a problem with %zu objectives",
                 inner_problem_->number_of_objectives());
  }

  if (inner_problem_->best_value().empty()) {
    /* Setting the value to 0 results in the assertion y>=optimal_fvalue being susceptible to failure */
    coco_warning("undefined optimal f value. Set to 0");
    optimal_fvalue_ = 0;
  } else {
    optimal_fvalue_ = inner_problem_->best_value()[0];
  }

  number_of_evaluations_ = 0;
  number_of_evaluations_constraints_ = 0;
  best_solution_.resize(inner_problem_->number_of_variables());
  /* TODO: the following inits are just to be in the safe side and
   * should eventually be removed. Some fields of the bbob_old_logger struct
   * might be useless
   */
  written_last_eval_ = 0;
  last_fvalue_ = DBL_MAX;
  is_initialized_ = 0;
  log_discrete_as_int_ = observer.log_discrete_as_int;

  /* Initialize triggers based on target values and number of evaluations */
  targets_ = ObserverLogTargets(observer.number_target_triggers, observer.log_target_precision);
  evaluation_triggers_ = ObserverEvaluations(observer.number_evaluation_triggers, observer.base_evaluation_triggers,
                                             inner_problem_->number_of_variables());
}

ProblemPtr ObserverBbobOld::AddLogger(ProblemPtr problem) {
  if (number_of_loggers() > 0)
    coco_error("The current bbob_old_logger (observer) must be closed before a new one is opened");
  return std::make_unique<LoggerBbobOld>(*this, std::move(problem));
}

} // namespace coco
