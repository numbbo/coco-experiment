#pragma once
#include "coco.h"
#include "coco_internal.h"

/**
 * @brief Data type for univariate function tasy_uv.
 */
typedef struct tasy_data_s {
  double beta;
  size_t i;
  size_t n;
} tasy_data;

/**
 * @brief Data type for transform_vars_asymmetric.
 */
typedef struct transform_vars_asymmetric_data_s {
  double *x;
  double beta;
} transform_vars_asymmetric_data_t;

void transform_inv_initial_asymmetric(coco_problem_t *problem, const double *xopt);
coco_problem_t *transform_vars_asymmetric(coco_problem_t *inner_problem, const double beta);
void transform_vars_asymmetric_free(void *thing);
void transform_vars_asymmetric_evaluate_constraint(coco_problem_t *problem, const double *x, double *y,
                                                   int update_counter);
void transform_vars_asymmetric_evaluate_function(coco_problem_t *problem, const double *x, double *y);
double tasy_uv(double xi, tasy_data *d);
double tasy_uv_inv(double yi, tasy_data *d);
