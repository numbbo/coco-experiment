#pragma once
#include "coco.h"
#include "coco_internal.h"
#undef INTERFACE
coco_problem_t *transform_vars_discretize(coco_problem_t *inner_problem, const double *smallest_values_of_interest,
                                          const double *largest_values_of_interest,
                                          const size_t number_of_integer_variables);
void transform_vars_discretize_free(void *thing);
void transform_vars_discretize_evaluate_function(coco_problem_t *problem, const double *x, double *y);
