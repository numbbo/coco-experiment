#pragma once
#include "coco.h"

coco_problem_t *f_linear_slope_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                     const size_t instance, const long rseed,
                                                     const char *problem_id_template,
                                                     const char *problem_name_template);
coco_problem_t *f_linear_slope_allocate(const size_t number_of_variables, const double *best_parameter);
void f_linear_slope_evaluate_gradient(coco_problem_t *problem, const double *x, double *y);
void f_linear_slope_evaluate(coco_problem_t *problem, const double *x, double *y);
double f_linear_slope_raw(const double *x, const size_t number_of_variables, const double *best_parameter);
