#pragma once
#include "coco.h"
#include "coco_observer.h"

#include <stdio.h>

/**
 * @brief The toy logger data type.
 */
typedef struct {
  FILE *log_file;                       /**< @brief Pointer to the file already prepared for logging. */
  coco_observer_log_targets_t *targets; /**< @brief Triggers based on logarithmic target values. */
  size_t number_of_evaluations;         /**< @brief The number of evaluations performed so far. */
  int precision_x;                      /**< @brief Precision for outputting decision values. */
  int precision_f;                      /**< @brief Precision for outputting objective values. */
} logger_toy_data_t;

coco_problem_t *logger_toy(coco_observer_t *observer, coco_problem_t *inner_problem);
void logger_toy_evaluate(coco_problem_t *problem, const double *x, double *y);
void logger_toy_free(void *stuff);
