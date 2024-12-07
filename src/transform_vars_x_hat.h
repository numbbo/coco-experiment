#pragma once
#include "coco.h"

coco_problem_t *transform_vars_x_hat(coco_problem_t *inner_problem, const long seed);
void transform_vars_x_hat_free(void *thing);
void transform_vars_x_hat_evaluate(coco_problem_t *problem, const double *x, double *y);
