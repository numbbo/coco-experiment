#pragma once
#include "coco.h"
#include "coco_internal.h"
#undef INTERFACE

/**
 * @brief Data type for the versatile_data_t
 */
typedef struct {
  coco_problem_t *sub_problem_mu0;
  coco_problem_t *sub_problem_mu1;
  double *x_hat;
} f_lunacek_bi_rastrigin_versatile_data_t;

coco_problem_t *transform_vars_x_hat_generic(coco_problem_t *inner_problem, const double *sign_vector);
void transform_vars_x_hat_generic_free(void *thing);
void transform_vars_x_hat_generic_evaluate(coco_problem_t *problem, const double *x, double *y);
