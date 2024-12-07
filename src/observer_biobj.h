#pragma once
#include "coco.h"
#include "coco_internal.h"

/** @brief Enum for denoting the way in which the nondominated solutions are treated. */
typedef enum { LOG_NONDOM_NONE, LOG_NONDOM_FINAL, LOG_NONDOM_ALL, LOG_NONDOM_READ } observer_biobj_log_nondom_e;

/** @brief Enum for denoting when the decision variables are logged. */
typedef enum { LOG_VARS_NEVER, LOG_VARS_LOW_DIM, LOG_VARS_ALWAYS } observer_biobj_log_vars_e;

/**
 * @brief The bbob-biobj observer data type.
 *
 * There is a cyclic reference between the biobj_logger and the biobj_observer (through the observer's
 * observed_problem and the logger's data, which points to the observer). This is needed to be able
 * to free both objects without problems.
 */
typedef struct observer_biobj_data_s {
  coco_problem_t *observed_problem; /**< @brief Pointer to the observed problem (NULL if none is observed) */

  observer_biobj_log_nondom_e log_nondom_mode; /**< @brief Handling of the nondominated solutions. */
  observer_biobj_log_vars_e log_vars_mode;     /**< @brief When the decision variables are logged. */

  int compute_indicators; /**< @brief Whether to compute indicators. */
  int produce_all_data;   /**< @brief Whether to produce all data. */
  int large_mdat_file_warning;    /**< @brief Whether the large mdat file warning has already been output. */

  long previous_function;  /**< @brief Function of the previous logged problem. */
  long previous_dimension; /**< @brief Dimension of the previous logged problem */

} observer_biobj_data_t;

void observer_biobj(coco_observer_t *observer, const char *options, coco_option_keys_t **option_keys);
void observer_biobj_data_free(void *stuff);
void logger_biobj_data_nullify_observer(void *logger_data);
void logger_biobj_signal_restart(coco_problem_t *problem);
void logger_biobj_free(void *logger);
coco_problem_t *logger_biobj(coco_observer_t *observer, coco_problem_t *problem);
