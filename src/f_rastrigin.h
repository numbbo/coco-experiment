#pragma once
#include "coco.h"
#undef INTERFACE
coco_problem_t *f_rastrigin_rotated_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                               const size_t instance, const long rseed,
                                                               const char *problem_id_template,
                                                               const char *problem_name_template);
coco_problem_t *f_rastrigin_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                       const size_t instance, const long rseed,
                                                       const char *problem_id_template,
                                                       const char *problem_name_template);
void f_rastrigin_cons_compute_xopt(double *xopt, const long rseed, const size_t dim);
coco_problem_t *f_rastrigin_permblockdiag_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                const size_t instance, const long rseed,
                                                                const char *problem_id_template,
                                                                const char *problem_name_template);
coco_problem_t *f_rastrigin_rotated_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                          const size_t instance, const long rseed,
                                                          const char *problem_id_template,
                                                          const char *problem_name_template);
coco_problem_t *f_rastrigin_bbob_problem_allocate(const size_t function, const size_t dimension, const size_t instance,
                                                  const long rseed, const char *problem_id_template,
                                                  const char *problem_name_template);
coco_problem_t *f_rastrigin_allocate(const size_t number_of_variables);
void f_rastrigin_evaluate_gradient(coco_problem_t *problem, const double *x, double *y);
void f_rastrigin_evaluate(coco_problem_t *problem, const double *x, double *y);
double f_rastrigin_raw(const double *x, const size_t number_of_variables);
