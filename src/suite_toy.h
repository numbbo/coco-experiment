#pragma once
#include "coco.h"
#undef INTERFACE
coco_problem_t *suite_toy_get_problem(coco_suite_t *suite, const size_t function_idx, const size_t dimension_idx,
                                      const size_t instance_idx);
coco_suite_t *suite_toy_initialize(void);
coco_suite_t *coco_suite_allocate(const char *suite_name, const size_t number_of_functions,
                                  const size_t number_of_dimensions, const size_t *dimensions,
                                  const char *default_instances, const int known_optima);
