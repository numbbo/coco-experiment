#pragma once
#include "coco.h"
#include "coco_internal.h"

/** @brief Enum for denoting when the decision variables and constraints are logged. */
typedef enum { LOG_NEVER, LOG_LOW_DIM, LOG_ALWAYS } observer_rw_log_e;

/**
 * @brief The real-world observer data type.
 */
typedef struct {
  observer_rw_log_e log_vars_mode; /**< @brief When the decision variables are logged. */
  observer_rw_log_e log_cons_mode; /**< @brief When the constraints are logged. */
  size_t low_dim_vars;             /**< @brief "Low dimension" for decision variables. */
  size_t low_dim_cons;             /**< @brief "Low dimension" for constraints. */
  int log_only_better;             /**< @brief Whether to log only solutions that are better than previous
                                               ones (only for the single-objective problems). */
  int log_time;                    /**< @brief Whether to log time. */
} observer_rw_data_t;

void observer_rw(coco_observer_t *observer, const char *options, coco_option_keys_t **option_keys);
void logger_rw_free(void *logger);
coco_problem_t *logger_rw(coco_observer_t *observer, coco_problem_t *problem);
