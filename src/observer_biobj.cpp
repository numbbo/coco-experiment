/**
 * @file observer_biobj.cpp
 * @brief Implementation of the bbob-biobj observer.
 */

#include "coco.h"
#include "coco_internal.h"

#include "coco_utilities.cpp"
#include "mo_utilities.cpp"

/** @brief Enum for denoting the way in which the nondominated solutions are treated. */
enum observer_biobj_log_nondom_e { LOG_NONDOM_NONE, LOG_NONDOM_FINAL, LOG_NONDOM_ALL, LOG_NONDOM_READ };

/** @brief Enum for denoting when the decision variables are logged. */
enum observer_biobj_log_vars_e { LOG_VARS_NEVER, LOG_VARS_LOW_DIM, LOG_VARS_ALWAYS };

struct logger_biobj_t;

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
struct observer_biobj_t final : public coco_observer_s {

  observer_biobj_t(const coco_observer_settings_s &observer_settings, std::string_view options,
                   coco_option_keys_t &option_keys)
      : coco_observer_s(observer_settings), observed_problem(nullptr), log_nondom_mode(LOG_NONDOM_ALL),
        log_vars_mode(LOG_VARS_LOW_DIM), compute_indicators(1), produce_all_data(0), large_mdat_file_warning(0),
        previous_function(-1), previous_dimension(-1) {

    std::string string_value;

    /* Sets the valid keys for bbob-biobj observer options
     * IMPORTANT: This list should be up-to-date with the code and the documentation */
    option_keys = {"log_nondominated", "log_decision_variables", "compute_indicators", "produce_all_data"};

    if (coco_options_read_string(options, "log_nondominated", string_value) > 0) {
      if (string_value == "none")
        log_nondom_mode = LOG_NONDOM_NONE;
      else if (string_value == "final")
        log_nondom_mode = LOG_NONDOM_FINAL;
      else if (string_value == "all")
        log_nondom_mode = LOG_NONDOM_ALL;
      else if (string_value == "read")
        log_nondom_mode = LOG_NONDOM_READ;
    }

    if (coco_options_read_string(options, "log_decision_variables", string_value) > 0) {
      if (string_value == "none")
        log_vars_mode = LOG_VARS_NEVER;
      else if (string_value == "all")
        log_vars_mode = LOG_VARS_ALWAYS;
      else if (string_value == "low_dim")
        log_vars_mode = LOG_VARS_LOW_DIM;
    }

    if (coco_options_read_int(options, "compute_indicators", &compute_indicators) == 0)
      compute_indicators = 1;

    if (coco_options_read_int(options, "produce_all_data", &produce_all_data) == 0)
      produce_all_data = 0;

    if (produce_all_data) {
      compute_indicators = 1;
      log_nondom_mode = LOG_NONDOM_ALL;
      log_vars_mode = LOG_VARS_LOW_DIM;
    }

    if ((log_nondom_mode == LOG_NONDOM_NONE) && (!compute_indicators)) {
      /* No logging required */
      is_active = 0;
    }
  }

  ~observer_biobj_t() override;

  coco_problem_ptr_t add_logger(coco_problem_ptr_t problem) override;

  void signal_restart(coco_problem_t *problem) override;

  logger_biobj_t *observed_problem; /**< @brief The logger of the observed problem (NULL if none is observed) */

  observer_biobj_log_nondom_e log_nondom_mode; /**< @brief Handling of the nondominated solutions. */
  observer_biobj_log_vars_e log_vars_mode;     /**< @brief When the decision variables are logged. */

  int compute_indicators;      /**< @brief Whether to compute indicators. */
  int produce_all_data;        /**< @brief Whether to produce all data. */
  int large_mdat_file_warning; /**< @brief Whether the large mdat file warning has already been output. */

  long previous_function;  /**< @brief Function of the previous logged problem. */
  long previous_dimension; /**< @brief Dimension of the previous logged problem */
};

