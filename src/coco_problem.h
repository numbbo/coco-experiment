#pragma once
#include "coco.h"
#include "coco_internal.h"

coco_problem_t *coco_problem_stacked_allocate(coco_problem_t *problem1, coco_problem_t *problem2,
                                              const double *smallest_values_of_interest,
                                              const double *largest_values_of_interest);
void coco_problem_stacked_free(coco_problem_t *problem);
void coco_problem_stacked_evaluate_constraint(coco_problem_t *problem, const double *x, double *y, int update_counter);
void coco_problem_stacked_evaluate_function(coco_problem_t *problem, const double *x, double *y);
coco_problem_t *coco_problem_transformed_allocate(coco_problem_t *inner_problem, void *user_data,
                                                  coco_data_free_function_t data_free_function,
                                                  const char *name_prefix);
void coco_problem_transformed_free(coco_problem_t *problem);
void coco_problem_transformed_free_data(coco_problem_t *problem);
void coco_problem_transformed_recommend_solution(coco_problem_t *problem, const double *x);
void bbob_problem_transformed_evaluate_gradient(coco_problem_t *problem, const double *x, double *y);
void coco_problem_transformed_evaluate_constraint(coco_problem_t *problem, const double *x, double *y,
                                                  int update_counter);
void coco_problem_transformed_evaluate_function(coco_problem_t *problem, const double *x, double *y);
coco_problem_t *coco_problem_transformed_get_inner_problem(const coco_problem_t *problem);
void *coco_problem_transformed_get_data(const coco_problem_t *problem);
void bbob_biobj_problem_best_parameter_print(const coco_problem_t *problem);
void logger_biobj_evaluate(coco_problem_t *problem, const double *x, double *y);
void bbob_problem_best_parameter_print(const coco_problem_t *problem);
int coco_problem_is_tainted(const coco_problem_t *problem);
int coco_problem_get_best_parameter(coco_problem_t *problem, double *best_parameter);
size_t coco_problem_get_suite_dep_instance(const coco_problem_t *problem);
size_t coco_problem_get_suite_dep_function(const coco_problem_t *problem);
size_t coco_problem_get_suite_dep_index(const coco_problem_t *problem);
void coco_problem_set_suite(coco_problem_t *problem, coco_suite_t *suite);
coco_suite_t *coco_problem_get_suite(const coco_problem_t *problem);
void coco_problem_get_initial_solution(const coco_problem_t *problem, double *initial_solution);
const double *coco_problem_get_largest_fvalues_of_interest(const coco_problem_t *problem);
size_t coco_problem_get_number_of_integer_variables(const coco_problem_t *problem);
const double *coco_problem_get_largest_values_of_interest(const coco_problem_t *problem);
const double *coco_problem_get_smallest_values_of_interest(const coco_problem_t *problem);
const char *coco_problem_get_type(const coco_problem_t *problem);
const char *coco_problem_get_id(const coco_problem_t *problem);
const char *coco_problem_get_name(const coco_problem_t *problem);
double coco_problem_get_final_target_fvalue1(const coco_problem_t *problem);
double coco_problem_get_best_value(const coco_problem_t *problem);
double coco_problem_get_best_observed_fvalue1(const coco_problem_t *problem);
int coco_problem_final_target_hit(const coco_problem_t *problem);
int coco_problem_best_parameter_not_zero(const coco_problem_t *problem);
size_t coco_problem_get_evaluations_constraints(const coco_problem_t *problem);
size_t coco_problem_get_evaluations(const coco_problem_t *problem);
void coco_problem_set_type(coco_problem_t *problem, const char *type, ...);
void coco_problem_set_name(coco_problem_t *problem, const char *name, ...);
void coco_problem_set_id(coco_problem_t *problem, const char *id, ...);
int coco_problem_id_is_fine(const char *id, ...);
void coco_problem_free(coco_problem_t *problem);
coco_problem_t *coco_problem_allocate_from_scalars(const char *problem_name, coco_evaluate_function_t evaluate_function,
                                                   coco_problem_free_function_t problem_free_function,
                                                   const size_t number_of_variables,
                                                   const double smallest_value_of_interest,
                                                   const double largest_value_of_interest, const double best_parameter);
coco_problem_t *coco_problem_duplicate(const coco_problem_t *other);
coco_problem_t *coco_problem_allocate(const size_t number_of_variables, const size_t number_of_objectives,
                                      const size_t number_of_constraints);
void coco_recommend_solution(coco_problem_t *problem, const double *x);
void bbob_evaluate_gradient(coco_problem_t *problem, const double *x, double *y);
void coco_evaluate_constraint(coco_problem_t *problem, const double *x, double *y);
void coco_evaluate_constraint_optional_update(coco_problem_t *problem, const double *x, double *y, int update_counter);
size_t coco_problem_get_number_of_constraints(const coco_problem_t *problem);
size_t coco_problem_get_number_of_objectives(const coco_problem_t *problem);
size_t coco_problem_get_dimension(const coco_problem_t *problem);
void coco_evaluate_function(coco_problem_t *problem, const double *x, double *y);
