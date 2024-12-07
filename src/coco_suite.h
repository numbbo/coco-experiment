#pragma once
#include "coco.h"


size_t coco_suite_encode_problem_index(const coco_suite_t *suite, const size_t function_idx, const size_t dimension_idx,
                                       const size_t instance_idx);
coco_problem_t *coco_suite_get_next_problem(coco_suite_t *suite, coco_observer_t *observer);
coco_suite_t *coco_suite(const char *suite_name, const char *suite_instance, const char *suite_options);
int coco_suite_is_next_dimension_found(coco_suite_t *suite);
int coco_suite_is_next_function_found(coco_suite_t *suite);
int coco_suite_is_next_instance_found(coco_suite_t *suite);
int coco_suite_is_next_item_found(const size_t *items, const size_t number_of_items, long *current_item_id);
size_t *coco_suite_get_instance_indices(const coco_suite_t *suite, const char *suite_instance);
size_t coco_suite_get_number_of_problems(const coco_suite_t *suite);
void coco_suite_decode_problem_index(const coco_suite_t *suite, const size_t problem_index, size_t *function_idx,
                                     size_t *dimension_idx, size_t *instance_idx);
coco_problem_t *coco_suite_get_problem(coco_suite_t *suite, const size_t problem_index);
void coco_suite_free(coco_suite_t *suite);
size_t coco_suite_get_instance_from_instance_index(const coco_suite_t *suite, const size_t instance_idx);
size_t coco_suite_get_dimension_from_dimension_index(const coco_suite_t *suite, const size_t dimension_idx);
size_t coco_suite_get_function_from_function_index(const coco_suite_t *suite, const size_t function_idx);
void coco_suite_filter_dimensions(coco_suite_t *suite, const size_t *dimension_numbers);
void coco_suite_filter_indices(size_t *items, const size_t number_of_items, const size_t *indices, const char *name);
void coco_suite_set_instance(coco_suite_t *suite, const size_t *instance_numbers);
coco_suite_t *coco_suite_allocate(const char *suite_name, const size_t number_of_functions,
                                  const size_t number_of_dimensions, const size_t *dimensions,
                                  const char *default_instances, const int known_optima);
coco_problem_t *coco_suite_get_problem_by_function_dimension_instance(coco_suite_t *suite, const size_t function,
                                                                      const size_t dimension, const size_t instance);
void coco_suite_get_best_indicator_value(const int known_optima, const coco_problem_t *problem,
                                         const char *indicator_name, double *value);
coco_problem_t *coco_suite_get_problem_from_indices(coco_suite_t *suite, const size_t function_idx,
                                                    const size_t dimension_idx, const size_t instance_idx);
const char *coco_suite_get_instances_by_year(const coco_suite_t *suite, const int year);
coco_suite_t *coco_suite_intialize(const char *suite_name);
