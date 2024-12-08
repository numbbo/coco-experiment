#pragma once
#include "coco.h"
#undef INTERFACE
coco_problem_t *transform_vars_shift(coco_problem_t *inner_problem, const double *offset,
                                     const int shift_constraint_only);
void transform_vars_shift_free(void *thing);
void transform_vars_shift_evaluate_gradient(coco_problem_t *problem, const double *x, double *y);
void transform_vars_shift_evaluate_constraint(coco_problem_t *problem, const double *x, double *y, int update_counter);
void transform_vars_shift_evaluate_function(coco_problem_t *problem, const double *x, double *y);
