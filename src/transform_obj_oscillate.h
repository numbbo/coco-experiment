#pragma once
#include "coco.h"

coco_problem_t *transform_obj_oscillate(coco_problem_t *inner_problem);
void transform_obj_oscillate_evaluate(coco_problem_t *problem, const double *x, double *y);
