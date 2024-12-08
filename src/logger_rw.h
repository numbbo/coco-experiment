#pragma once
#include "coco.h"

#include <stdio.h>

/**
 * @brief The rw logger data type.
 *
 * @note Some fields from the observers (coco_observer as well as observer_rw) need to be copied here
 * because the observers can be deleted before the logger is finalized and we need these fields for
 * finalization.
 */
typedef struct {
  FILE *out_file;              /**< @brief File for logging. */
  size_t num_func_evaluations; /**< @brief The number of function evaluations performed so far. */
  size_t num_cons_evaluations; /**< @brief The number of evaluations of constraints performed so far. */

  double best_value;    /**< @brief The best-so-far value. */
  double current_value; /**< @brief The current value. */

  int log_vars;        /**< @brief Whether to log the decision values. */
  int log_cons;        /**< @brief Whether to log the constraints. */
  int log_only_better; /**< @brief Whether to log only solutions that are better than previous ones. */
  int log_time;        /**< @brief Whether to log evaluation time. */

  int precision_x;         /**< @brief Precision for outputting decision values. */
  int precision_f;         /**< @brief Precision for outputting objective values. */
  int precision_g;         /**< @brief Precision for outputting constraint values. */
  int log_discrete_as_int; /**< @brief Whether to output discrete variables in int or double format. */
} logger_rw_data_t;

coco_problem_t *logger_rw(coco_observer_t *observer, coco_problem_t *inner_problem);
void logger_rw_free(void *stuff);
void logger_rw_evaluate(coco_problem_t *problem, const double *x, double *y);
