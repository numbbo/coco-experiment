#pragma once
#include "coco.h"
#include "coco_internal.h"
#undef INTERFACE
coco_problem_t *suite_bbob_mixint_get_problem(coco_suite_t *suite, const size_t function_idx,
                                              const size_t dimension_idx, const size_t instance_idx);
coco_problem_t *coco_get_bbob_mixint_problem(const size_t function, const size_t dimension, const size_t instance,
                                             const coco_get_problem_function_t coco_get_problem_function);
const char *suite_bbob_mixint_get_instances_by_year(const int year);
coco_suite_t *suite_bbob_mixint_initialize(const char *suite_name);
extern double suite_bbob_mixint_scaling_factors[];
coco_suite_t *coco_suite_allocate(const char *suite_name, const size_t number_of_functions,
                                  const size_t number_of_dimensions, const size_t *dimensions,
                                  const char *default_instances, const int known_optima);
