#pragma once
#include "coco.h"
#undef INTERFACE
coco_problem_t *transform_vars_enforce_box(coco_problem_t *inner_problem);
void transform_vars_enforce_box_evaluate_function(coco_problem_t *problem, const double *x, double *y);
