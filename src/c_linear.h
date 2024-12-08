#pragma once
#include "coco.h"
#undef INTERFACE
void con_update_linear_combination(double *linear_combination, const coco_problem_t *problem, double weight);
coco_problem_t *c_linear_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                    const size_t instance, const size_t number_of_active_constraints,
                                                    const char *problem_id_template, const char *problem_name_template,
                                                    const double *feasible_direction);
coco_problem_t *c_linear_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                    const size_t instance, const size_t number_of_active_constraints,
                                                    const char *problem_id_template, const char *problem_name_template,
                                                    const double *feasible_direction);
coco_problem_t *c_linear_single_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                           const size_t instance, const size_t constraint_number,
                                                           const double factor1, const char *problem_id_template,
                                                           const char *problem_name_template, double *gradient,
                                                           double x_shift_factor, const double *feasible_direction);
coco_problem_t *c_linear_single_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                           const size_t instance, const size_t constraint_number,
                                                           const double factor1, const char *problem_id_template,
                                                           const char *problem_name_template, double *gradient,
                                                           const double x_shift_factor,
                                                           const double *feasible_direction);
coco_problem_t *c_linear_transform(coco_problem_t *inner_problem, const double *gradient, const double x_shift_factor);
coco_problem_t *c_linear_transform(coco_problem_t *inner_problem, const double *gradient, const double shift_factor);
coco_problem_t *c_sum_variables_allocate(const size_t number_of_variables);
coco_problem_t *c_sum_variables_allocate(const size_t number_of_variables);
void c_linear_gradient_free(void *thing);
void c_linear_gradient_free(void *thing);
coco_problem_t *c_guarantee_feasible_point(coco_problem_t *problem, const double *feasible_point);
coco_problem_t *c_guarantee_feasible_point(coco_problem_t *problem, const double *feasible_direction);
void c_linear_single_evaluate(coco_problem_t *self, const double *x, double *y, int update_counter);
void c_linear_single_evaluate(coco_problem_t *self, const double *x, double *y, int update_counter);
void c_sum_variables_evaluate(coco_problem_t *self, const double *x, double *y, int update_counter);
void c_sum_variables_evaluate(coco_problem_t *self, const double *x, double *y, int update_counter);
