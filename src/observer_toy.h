#pragma once
#include "coco.h"
#include "coco_internal.h"

#include <stdio.h>

/**
 * @brief The toy observer data type.
 */
typedef struct {
  FILE *log_file; /**< @brief File used for logging. */
} observer_toy_data_t;

void observer_toy(coco_observer_t *observer, const char *options, coco_option_keys_t **option_keys);
void observer_toy_free(void *stuff);
void logger_toy_free(void *logger);
coco_problem_t *logger_toy(coco_observer_t *observer, coco_problem_t *problem);
