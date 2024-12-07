#pragma once
#include "coco.h"

/**
 * @brief Data type in problem->versatile_data of f_gallagher.c
 */
typedef struct f_gallagher_versatile_data_s {
  double *x;
  size_t number_of_peaks;
  coco_problem_t **sub_problems;
  double *rotated_x;
  size_t nb_blocks, *block_sizes, *block_size_map, *first_non_zero_map;
  double **B;
} f_gallagher_versatile_data_t;

/**
 * @brief Data type for transform_vars_gallagher_blockrotation.
 */
typedef struct transform_vars_gallagher_blockrotation_t{
  double *x;
} transform_vars_gallagher_blockrotation_t;

COCO_UNUSED coco_problem_t *transform_vars_gallagher_blockrotation(coco_problem_t *inner_problem);
coco_problem_t *transform_vars_gallagher_blockrotation(coco_problem_t *inner_problem);
void transform_vars_gallagher_blockrotation_evaluate(coco_problem_t *problem, const double *x, double *y);
void f_gallagher_versatile_data_free(coco_problem_t *problem);
void transform_vars_gallagher_blockrotation_free(void *thing);
