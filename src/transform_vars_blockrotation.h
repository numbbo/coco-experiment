#pragma once
#include "coco.h"

coco_problem_t *transform_vars_blockrotation(coco_problem_t *inner_problem, const double *const *B,
                                             const size_t number_of_variables, const size_t *block_sizes,
                                             const size_t nb_blocks);
void transform_vars_blockrotation_test(coco_problem_t *problem, double precision);
void transform_vars_blockrotation_free(void *stuff);
void transform_vars_blockrotation_evaluate(coco_problem_t *problem, const double *x, double *y);
void transform_vars_blockrotation_apply(coco_problem_t *problem, const double *x, double *y);
void transform_vars_blockrotation_get_row(coco_problem_t *problem, size_t i, double *y);
