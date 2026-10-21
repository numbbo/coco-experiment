/**
 * @file observer_biobj.cpp
 * @brief Implementation of the bbob-biobj observer.
 */

#include "coco.h"
#include "coco_internal.h"

#include "coco_utilities.cpp"
#include "mo_utilities.cpp"

namespace coco {

/** @brief Enum for denoting the way in which the nondominated solutions are treated. */
enum ObserverBiobjLogNondom { LOG_NONDOM_NONE, LOG_NONDOM_FINAL, LOG_NONDOM_ALL, LOG_NONDOM_READ };

/** @brief Enum for denoting when the decision variables are logged. */
enum ObserverBiobjLogVars { LOG_VARS_NEVER, LOG_VARS_LOW_DIM, LOG_VARS_ALWAYS };

class LoggerBiobj;

/**
 * @brief The bi-objective observer.
 *
 * There is a cyclic reference between the biobj_logger and the biobj_observer (through the observer's
 * observed_problem and the logger's data, which points to the observer). This is needed to be able
 * to free both objects without problems.
 *
 * Possible options:
 *
 * - "log_nondominated: STRING" determines how the nondominated solutions are handled. STRING can take on the
 * values "none" (don't log nondominated solutions), "final" (log only the final nondominated solutions),
 * "all" (log every solution that is nondominated at creation time) and "read" (the nondominated solutions
 * are not logged, but are passed to the logger as input - this is a functionality needed in pre-processing
 * of the data). The default value is "all".
 *
 * - "log_decision_variables: STRING" determines whether the decision variables are to be logged in addition
 * to the objective variables in the output of nondominated solutions. STRING can take on the values "none"
 * (don't output decision variables), "low_dim"(output decision variables only for dimensions lower or equal
 * to 5) and "all" (output all decision variables). The default value is "low_dim".
 *
 * - "compute_indicators: VALUE" determines whether to compute and output performance indicators (1) or not
 * (0). The default value is 1.
 *
 * - "produce_all_data: VALUE" determines whether to produce all data required for the workshop. If set to 1,
 * it overwrites some other options and is equivalent to setting "log_nondominated: all",
 * "log_decision_variables: low_dim" and "compute_indicators: 1". If set to 0, it does not change the values
 * of the other options. The default value is 0.
 */
class ObserverBiobj final : public Observer {
public:
  ObserverBiobj(const ObserverSettings &observer_settings, std::string_view options, OptionKeys &option_keys)
      : Observer(observer_settings), observed_problem_(nullptr), log_nondom_mode_(LOG_NONDOM_ALL),
        log_vars_mode_(LOG_VARS_LOW_DIM), compute_indicators_(true), produce_all_data_(false),
        large_mdat_file_warning_(false), previous_function_(-1), previous_dimension_(-1) {

    std::string string_value;

    /* Sets the valid keys for bbob-biobj observer options
     * IMPORTANT: This list should be up-to-date with the code and the documentation */
    option_keys = {"log_nondominated", "log_decision_variables", "compute_indicators", "produce_all_data"};

    if (OptionsReadString(options, "log_nondominated", string_value) > 0) {
      if (string_value == "none")
        log_nondom_mode_ = LOG_NONDOM_NONE;
      else if (string_value == "final")
        log_nondom_mode_ = LOG_NONDOM_FINAL;
      else if (string_value == "all")
        log_nondom_mode_ = LOG_NONDOM_ALL;
      else if (string_value == "read")
        log_nondom_mode_ = LOG_NONDOM_READ;
    }

    if (OptionsReadString(options, "log_decision_variables", string_value) > 0) {
      if (string_value == "none")
        log_vars_mode_ = LOG_VARS_NEVER;
      else if (string_value == "all")
        log_vars_mode_ = LOG_VARS_ALWAYS;
      else if (string_value == "low_dim")
        log_vars_mode_ = LOG_VARS_LOW_DIM;
    }

    if (OptionsReadFlag(options, "compute_indicators", &compute_indicators_) == 0)
      compute_indicators_ = true;

    if (OptionsReadFlag(options, "produce_all_data", &produce_all_data_) == 0)
      produce_all_data_ = false;

    if (produce_all_data_) {
      compute_indicators_ = 1;
      log_nondom_mode_ = LOG_NONDOM_ALL;
      log_vars_mode_ = LOG_VARS_LOW_DIM;
    }

    if ((log_nondom_mode_ == LOG_NONDOM_NONE) && (!compute_indicators_)) {
      /* No logging required */
      Deactivate();
    }
  }

  ~ObserverBiobj() override;

  ProblemPtr AddLogger(ProblemPtr problem) override;

  void SignalRestart(Problem *problem) override;

private:
  /* The logger and its indicators read the settings and keep the state of the last logged problem */
  friend class LoggerBiobj;
  friend class LoggerBiobjIndicator;

  LoggerBiobj *observed_problem_; /**< @brief The logger of the observed problem (nullptr if none is observed) */
  ObserverBiobjLogNondom log_nondom_mode_;      /**< @brief Handling of the nondominated solutions. */
  ObserverBiobjLogVars log_vars_mode_;          /**< @brief When the decision variables are logged. */
  bool compute_indicators_;                     /**< @brief Whether to compute indicators. */
  bool produce_all_data_;                       /**< @brief Whether to produce all data. */
  bool large_mdat_file_warning_; /**< @brief Whether the large mdat file warning has already been output. */
  long previous_function_;       /**< @brief Function of the previous logged problem. */
  long previous_dimension_;      /**< @brief Dimension of the previous logged problem */
};

} // namespace coco
