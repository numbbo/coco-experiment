#pragma once
#include "coco.h"
#include "coco_internal.h"

coco_problem_t *transform_obj_gaussian_noise(coco_problem_t *inner_problem, const double beta);
void transform_obj_gaussian_noise_evaluate_function(coco_problem_t *problem, const double *x, double *y);
