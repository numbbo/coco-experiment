#pragma once
#include "coco.h"
#include "coco_internal.h"
#undef INTERFACE
coco_problem_t *transform_obj_penalize(coco_problem_t *inner_problem, const double factor);
void transform_obj_penalize_evaluate(coco_problem_t *problem, const double *x, double *y);
