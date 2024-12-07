#pragma once
#include "coco.h"

coco_problem_t *transform_vars_permblockdiag(coco_problem_t *inner_problem, const double *const *B, const size_t *P1,
                                             const size_t *P2, const size_t number_of_variables,
                                             const size_t *block_sizes, const size_t nb_blocks);
void transform_vars_permblockdiag_free(void *thing);
void transform_vars_permblockdiag_evaluate(coco_problem_t *problem, const double *x, double *y);
