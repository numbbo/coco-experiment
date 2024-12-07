#pragma once
#include "coco.h"

coco_problem_t *transform_vars_scale(coco_problem_t *inner_problem, const double factor);
void transform_vars_scale_free(void *thing);
void transform_vars_scale_evaluate(coco_problem_t *problem, const double *x, double *y);
