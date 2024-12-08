#pragma once
#include "coco.h"
#undef INTERFACE
coco_problem_t *transform_obj_uniform_noise(coco_problem_t *inner_problem, const double alpha, const double beta);
void transform_obj_uniform_noise_evaluate_function(coco_problem_t *problem, const double *x, double *y);
