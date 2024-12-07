#pragma once
#include "coco.h"

#include "transform_vars_permutation_helpers.h"

coco_problem_t *transform_vars_inverse_permutation(coco_problem_t *inner_problem, const size_t *P,
                                                   const size_t number_of_variables);
void transform_vars_inverse_permutation_evaluate(coco_problem_t *problem, const double *x, double *y);
coco_problem_t *transform_vars_permutation(coco_problem_t *inner_problem, const size_t *P,
                                           const size_t number_of_variables);
void transform_vars_permutation_free(void *thing);
void transform_vars_permutation_evaluate(coco_problem_t *problem, const double *x, double *y);
