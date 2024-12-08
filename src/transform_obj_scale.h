#pragma once
#include "coco.h"
#undef INTERFACE
coco_problem_t *transform_obj_scale(coco_problem_t *inner_problem, const double factor);
void transform_obj_scale_evaluate_gradient(coco_problem_t *problem, const double *x, double *y);
void transform_obj_scale_evaluate_function(coco_problem_t *problem, const double *x, double *y);
