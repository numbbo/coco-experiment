#pragma once
#include "coco.h"

typedef struct f_step_ellipsoid_data_s f_step_ellipsoid_data_t;
typedef struct f_step_ellipsoid_versatile_data_s f_step_ellipsoid_versatile_data_t;

coco_problem_t *f_step_ellipsoid_permblockdiag_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                     const size_t instance, const long rseed,
                                                                     const char *problem_id_template,
                                                                     const char *problem_name_template);
coco_problem_t *f_step_ellipsoid_allocate(const size_t number_of_variables);
void f_step_ellipsoid_versatile_data_free(coco_problem_t *problem);
void f_step_ellipsoid_permblock_evaluate(coco_problem_t *problem, const double *x, double *y);
double f_step_ellipsoid_core(const double *x, const size_t number_of_variables,
                             f_step_ellipsoid_versatile_data_t *f_step_ellipsoid_versatile_data);
coco_problem_t *f_step_ellipsoid_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                       const size_t instance, const long rseed, const void *args,
                                                       const char *problem_id_template,
                                                       const char *problem_name_template);
void f_step_ellipsoid_free(coco_problem_t *problem);
void f_step_ellipsoid_evaluate(coco_problem_t *problem, const double *x, double *y);
double f_step_ellipsoid_raw(const double *x, const size_t number_of_variables, const f_step_ellipsoid_data_t *data);
