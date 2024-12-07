#pragma once
#include "coco.h"
#include "coco_internal.h"

coco_problem_t *transform_obj_cauchy_noise(coco_problem_t *inner_problem, const double alpha, const double p);
void transform_obj_cauchy_noise_evaluate_function(coco_problem_t *problem, const double *x, double *y);
