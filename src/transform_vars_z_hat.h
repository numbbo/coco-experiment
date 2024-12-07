#pragma once
#include "coco.h"

coco_problem_t *transform_vars_z_hat(coco_problem_t *inner_problem, const double *xopt);
void transform_vars_z_hat_free(void *thing);
void transform_vars_z_hat_evaluate(coco_problem_t *problem, const double *x, double *y);
