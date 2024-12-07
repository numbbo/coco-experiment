#pragma once
#include "coco.h"

#include "transform_vars_gallagher_blockrotation.h"

typedef struct f_gallagher_data_s f_gallagher_data_t;
typedef struct f_gallagher_permutation_s f_gallagher_permutation_t;

coco_problem_t *f_gallagher_permblockdiag_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                const size_t instance, const long rseed,
                                                                size_t number_of_peaks, const char *problem_id_template,
                                                                const char *problem_name_template);
coco_problem_t *f_gallagher_problem_allocate(const size_t number_of_variables, size_t number_of_peaks);
void f_gallagher_evaluate_core(coco_problem_t *problem, const double *x, double *y);
double f_gallagher_core(const double *x, size_t number_of_variables, f_gallagher_versatile_data_t *versatile_data);
coco_problem_t *f_gallagher_sub_problem_allocate(const size_t number_of_variables);
void f_gallagher_sub_evaluate_core(coco_problem_t *problem_i, const double *x, double *y);
coco_problem_t *f_gallagher_bbob_problem_allocate(const size_t function, const size_t dimension, const size_t instance,
                                                  const long rseed, const void *args, const char *problem_id_template,
                                                  const char *problem_name_template);
void f_gallagher_free(coco_problem_t *problem);
void f_gallagher_evaluate(coco_problem_t *problem, const double *x, double *y);
double f_gallagher_raw(const double *x, const size_t number_of_variables, f_gallagher_data_t *data);
int f_gallagher_compare_doubles(const void *a, const void *b);
