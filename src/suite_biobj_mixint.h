#pragma once
#include "coco.h"
#include "coco_internal.h"

#include "suite_biobj_utilities.h"

// typedef struct suite_biobj_new_inst_s suite_biobj_new_inst_t;

coco_problem_t *suite_biobj_mixint_get_problem(coco_suite_t *suite, const size_t function_idx,
                                               const size_t dimension_idx, const size_t instance_idx);
coco_problem_t *coco_get_biobj_mixint_problem(const size_t function, const size_t dimension, const size_t instance,
                                              const coco_get_problem_function_t coco_get_problem_function,
                                              suite_biobj_new_inst_t **new_inst_data, const size_t num_new_instances,
                                              const size_t *dimensions, const size_t num_dimensions);
const char *suite_biobj_mixint_get_instances_by_year(const int year);
coco_suite_t *suite_biobj_mixint_initialize(void);
void suite_biobj_new_inst_free(void *stuff);
coco_suite_t *coco_suite_allocate(const char *suite_name, const size_t number_of_functions,
                                  const size_t number_of_dimensions, const size_t *dimensions,
                                  const char *default_instances, const int known_optima);
