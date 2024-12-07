#pragma once
#include "coco.h"

#include "transform_vars_x_hat_generic.h"

typedef struct f_lunacek_bi_rastrigin_data_s f_lunacek_bi_rastrigin_data_t;

coco_problem_t *f_lunacek_bi_rastrigin_permblockdiag_bbob_problem_allocate(const size_t function,
                                                                           const size_t dimension,
                                                                           const size_t instance, const long rseed,
                                                                           const char *problem_id_template,
                                                                           const char *problem_name_template);
coco_problem_t *f_lunacek_bi_rastrigin_problem_allocate(const size_t number_of_variables);
void f_lunacek_bi_rastrigin_evaluate_core(coco_problem_t *problem, const double *x, double *y);
double f_lunacek_bi_rastrigin_core(const double *x, const size_t number_of_variables,
                                   f_lunacek_bi_rastrigin_versatile_data_t *f_lunacek_bi_rastrigin_versatile_data);
coco_problem_t *f_lunacek_bi_rastrigin_sub_problem_allocate(const size_t number_of_variables);
void f_lunacek_bi_rastrigin_sub_evaluate_core(coco_problem_t *problem, const double *x, double *y);
void f_lunacek_bi_rastrigin_versatile_data_free(coco_problem_t *problem);
coco_problem_t *f_lunacek_bi_rastrigin_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                             const size_t instance, const long rseed,
                                                             const char *problem_id_template,
                                                             const char *problem_name_template);
void f_lunacek_bi_rastrigin_free(coco_problem_t *problem);
void f_lunacek_bi_rastrigin_evaluate(coco_problem_t *problem, const double *x, double *y);
double f_lunacek_bi_rastrigin_raw(const double *x, const size_t number_of_variables,
                                  f_lunacek_bi_rastrigin_data_t *data);
