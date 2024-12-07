#pragma once
#include "coco.h"
#include "coco_internal.h"
#undef INTERFACE
coco_problem_t *transform_vars_conditioning(coco_problem_t *inner_problem, const double alpha);
void transform_vars_conditioning_free(void *thing);
void transform_vars_conditioning_evaluate_gradient(coco_problem_t *problem, const double *x, double *y);
void transform_vars_conditioning_evaluate(coco_problem_t *problem, const double *x, double *y);
