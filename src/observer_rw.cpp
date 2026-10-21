/**
 * @file observer_rw.cpp
 * @brief Implementation of an observer for real-world problems.
 */

#include "coco.h"
#include "coco_internal.h"

#include "coco_utilities.cpp"

/** @brief Enum for denoting when the decision variables and constraints are logged. */
enum observer_rw_log_e { LOG_NEVER, LOG_LOW_DIM, LOG_ALWAYS };

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
struct observer_rw_t final : public coco_observer_s {

  observer_rw_t(const coco_observer_settings_s &observer_settings, std::string_view options,
                coco_option_keys_t &option_keys)
      : coco_observer_s(observer_settings), log_vars_mode(LOG_ALWAYS), log_cons_mode(LOG_ALWAYS) {

    std::string string_value;

    /* Sets the valid keys for rw observer options
     * IMPORTANT: This list should be up-to-date with the code and the documentation */
    option_keys = {"log_variables", "log_constraints", "low_dim_vars", "low_dim_cons", "log_only_better", "log_time"};

    if (coco_options_read_string(options, "log_variables", string_value) > 0) {
      if (string_value == "none")
        log_vars_mode = LOG_NEVER;
      else if (string_value == "all")
        log_vars_mode = LOG_ALWAYS;
      else if (string_value == "low_dim")
        log_vars_mode = LOG_LOW_DIM;
    }

    if (coco_options_read_string(options, "log_constraints", string_value) > 0) {
      if (string_value == "none")
        log_cons_mode = LOG_NEVER;
      else if (string_value == "all")
        log_cons_mode = LOG_ALWAYS;
      else if (string_value == "low_dim")
        log_cons_mode = LOG_LOW_DIM;
    }

    if (coco_options_read_size_t(options, "low_dim_vars", &low_dim_vars) == 0)
      low_dim_vars = 10;

    if (coco_options_read_size_t(options, "low_dim_cons", &low_dim_cons) == 0)
      low_dim_cons = 10;

    if (coco_options_read_flag(options, "log_only_better", &log_only_better) == 0)
      log_only_better = false;

    if (coco_options_read_flag(options, "log_time", &log_time) == 0)
      log_time = false;
  }

  coco_problem_ptr_t add_logger(coco_problem_ptr_t problem) override;

  observer_rw_log_e log_vars_mode; /**< @brief When the decision variables are logged. */
  observer_rw_log_e log_cons_mode; /**< @brief When the constraints are logged. */
  size_t low_dim_vars;             /**< @brief "Low dimension" for decision variables. */
  size_t low_dim_cons;             /**< @brief "Low dimension" for constraints. */
  bool log_only_better = false;    /**< @brief Whether to log only solutions that are better than previous
                                      ones (only for the single-objective problems). */
  bool log_time = false;           /**< @brief Whether to log time. */
};

