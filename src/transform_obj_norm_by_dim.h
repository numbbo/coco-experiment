#pragma once
#include "coco.h"
#include "coco_internal.h"

coco_problem_t *transform_obj_norm_by_dim(coco_problem_t *inner_problem);
void transform_obj_norm_by_dim_evaluate(coco_problem_t *problem, const double *x, double *y);
