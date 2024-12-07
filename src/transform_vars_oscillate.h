#pragma once
#include "coco.h"

/**
 * @brief Data type for univariate function tosz_uv
 */
typedef struct tosz_data_s {
  double alpha;
} tosz_data;

/**
 * @brief Data type for transform_vars_oscillate.
 */
typedef struct transform_vars_oscillate_data_s {
  double alpha;
  double *oscillated_x;
} transform_vars_oscillate_data_t;

void transform_inv_initial_oscillate(coco_problem_t *problem, const double *xopt);
coco_problem_t *transform_vars_oscillate(coco_problem_t *inner_problem);
void transform_vars_oscillate_free(void *thing);
void transform_vars_oscillate_evaluate_constraint(coco_problem_t *problem, const double *x, double *y,
                                                  int update_counter);
void transform_vars_oscillate_evaluate_function(coco_problem_t *problem, const double *x, double *y);
double tosz_uv(double xi, tosz_data *d);
double tosz_uv_inv(double yi, tosz_data *d);
