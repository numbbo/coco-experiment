#pragma once
#include "coco.h"
#include "coco_internal.h"
#undef INTERFACE
coco_problem_t *transform_obj_power(coco_problem_t *inner_problem, const double exponent);
void transform_obj_power_evaluate(coco_problem_t *problem, const double *x, double *y);
