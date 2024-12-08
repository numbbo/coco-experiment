#pragma once
#include "coco.h"

coco_problem_t *transform_vars_affine(coco_problem_t *inner_problem, const double *M, const double *b,
                                      const size_t number_of_variables);
void transform_vars_affine_free(void *thing);
void transform_vars_affine_evaluate_gradient(coco_problem_t *problem, const double *x, double *y);
void transform_vars_affine_evaluate_constraint(coco_problem_t *problem, const double *x, double *y, int update_counter);
void transform_vars_affine_evaluate_function(coco_problem_t *problem, const double *x, double *y);
