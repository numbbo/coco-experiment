/**
 * @file observer_rw.cpp
 * @brief Implementation of an observer for real-world problems.
 */

#include "coco.h"
#include "coco_internal.h"

#include "coco_utilities.cpp"

namespace coco {

/** @brief Enum for denoting when the decision variables and constraints are logged. */
enum ObserverRwLog { LOG_NEVER, LOG_LOW_DIM, LOG_ALWAYS };

/**
 * @brief The observer for real-world problems.
 *
 * Possible options:
 *
 * - "log_variables: STRING" determines whether the decision variables are to be logged. STRING can take on
 * the values "none" (don't output decision variables), "low_dim"(output decision variables only for
 * dimensions lower or equal to low_dim_vars) and "all" (output all decision variables). The default value
 * is "all".
 *
 * - "log_constraints: STRING" determines whether the constraints are to be logged. STRING can take on the
 * values "none" (don't output constraints), "low_dim"(output constraints only for dimensions lower or equal
 * to low_dim_cons) and "all" (output all constraints). The default value is "all".
 *
 * - "low_dim_vars: VALUE" determines the value used to define "low_dim" for decision variables. The default
 * value is 10.
 *
 * - "low_dim_cons: VALUE" determines the value used to define "low_dim" for constraints. The default value
 * is 10.
 *
 * - "log_only_better: 0/1" determines whether all solutions are logged (0, default) or only the ones that
 * are better than previous ones (1). This is applicable only for the single-objective problems, where the
 * default value is 0. For multi-objective problems, all solutions are always logged.
 *
 * - "log_time: 0/1" determines whether the time needed to evaluate each solution is logged (0) or not (1).
 * The default value is 0.
 */
class ObserverRw final : public Observer {
public:
  ObserverRw(const ObserverSettings &observer_settings, std::string_view options, OptionKeys &option_keys)
      : Observer(observer_settings), log_vars_mode_(LOG_ALWAYS), log_cons_mode_(LOG_ALWAYS) {

    std::string string_value;

    /* Sets the valid keys for rw observer options
     * IMPORTANT: This list should be up-to-date with the code and the documentation */
    option_keys = {"log_variables", "log_constraints", "low_dim_vars", "low_dim_cons", "log_only_better", "log_time"};

    if (OptionsReadString(options, "log_variables", string_value) > 0) {
      if (string_value == "none")
        log_vars_mode_ = LOG_NEVER;
      else if (string_value == "all")
        log_vars_mode_ = LOG_ALWAYS;
      else if (string_value == "low_dim")
        log_vars_mode_ = LOG_LOW_DIM;
    }

    if (OptionsReadString(options, "log_constraints", string_value) > 0) {
      if (string_value == "none")
        log_cons_mode_ = LOG_NEVER;
      else if (string_value == "all")
        log_cons_mode_ = LOG_ALWAYS;
      else if (string_value == "low_dim")
        log_cons_mode_ = LOG_LOW_DIM;
    }

    if (OptionsReadSizeT(options, "low_dim_vars", &low_dim_vars_) == 0)
      low_dim_vars_ = 10;

    if (OptionsReadSizeT(options, "low_dim_cons", &low_dim_cons_) == 0)
      low_dim_cons_ = 10;

    if (OptionsReadFlag(options, "log_only_better", &log_only_better_) == 0)
      log_only_better_ = false;

    if (OptionsReadFlag(options, "log_time", &log_time_) == 0)
      log_time_ = false;
  }

  ProblemPtr AddLogger(ProblemPtr problem) override;

private:
  friend class LoggerRw; /* Reads the settings */

  ObserverRwLog log_vars_mode_;     /**< @brief When the decision variables are logged. */
  ObserverRwLog log_cons_mode_;     /**< @brief When the constraints are logged. */
  size_t low_dim_vars_;             /**< @brief "Low dimension" for decision variables. */
  size_t low_dim_cons_;             /**< @brief "Low dimension" for constraints. */
  bool log_only_better_ = false;    /**< @brief Whether to log only solutions that are better than previous
                                      ones (only for the single-objective problems). */
  bool log_time_ = false;           /**< @brief Whether to log time. */
};

} // namespace coco
