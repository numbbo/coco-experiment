#pragma once
#include "coco.h"
#undef INTERFACE
coco_problem_t *f_schwefel_generalized_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                             const size_t instance, const long rseed,
                                                             const char *problem_id_template,
                                                             const char *problem_name_template);
coco_problem_t *f_schwefel_generalized_allocate(const size_t number_of_variables);
void f_schwefel_generalized_evaluate(coco_problem_t *problem, const double *x, double *y);
double f_schwefel_generalized_raw(const double *x, const size_t number_of_variables);
