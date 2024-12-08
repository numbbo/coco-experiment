#pragma once
#include "coco.h"

coco_problem_t *transform_obj_shift(coco_problem_t *inner_problem, const double offset);
void transform_obj_shift_evaluate_gradient(coco_problem_t *problem, const double *x, double *y);
void transform_obj_shift_evaluate_function(coco_problem_t *problem, const double *x, double *y);
