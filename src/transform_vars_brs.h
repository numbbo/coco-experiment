#pragma once
#include "coco.h"
#undef INTERFACE
coco_problem_t *transform_vars_brs(coco_problem_t *inner_problem);
void transform_vars_brs_free(void *thing);
void transform_vars_brs_evaluate(coco_problem_t *problem, const double *x, double *y);
