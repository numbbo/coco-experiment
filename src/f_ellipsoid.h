#pragma once
#include "coco.h"

typedef struct f_ellipsoid_data_s f_ellipsoid_data_t;

coco_problem_t *f_ellipsoid_rotated_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                               const size_t instance, const long rseed,
                                                               const char *problem_id_template,
                                                               const char *problem_name_template);
coco_problem_t *f_ellipsoid_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                       const size_t instance, const long rseed,
                                                       const char *problem_id_template,
                                                       const char *problem_name_template);
coco_problem_t *f_ellipsoid_permblockdiag_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                const size_t instance, const long rseed,
                                                                const char *problem_id_template,
                                                                const char *problem_name_template);
coco_problem_t *f_ellipsoid_rotated_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                          const size_t instance, const long rseed, const void *args,
                                                          const char *problem_id_template,
                                                          const char *problem_name_template);
coco_problem_t *f_ellipsoid_bbob_problem_allocate(const size_t function, const size_t dimension, const size_t instance,
                                                  const long rseed, const char *problem_id_template,
                                                  const char *problem_name_template);
coco_problem_t *f_ellipsoid_allocate(const size_t number_of_variables, const double conditioning);
void f_ellipsoid_evaluate_gradient(coco_problem_t *problem, const double *x, double *y);
void f_ellipsoid_evaluate(coco_problem_t *problem, const double *x, double *y);
double f_ellipsoid_raw(const double *x, const size_t number_of_variables, f_ellipsoid_data_t *data);
