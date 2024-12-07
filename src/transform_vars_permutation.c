/**
 * @file transform_vars_permutation.c
 * @brief Implementation of permuting the decision values.
 */
#include "transform_vars_permutation.h"

#include <assert.h>

#include "coco.h"
#include "coco_utilities.h"
#include "transform_vars_permutation_helpers.h"

/**
 * @brief Data type for transform_vars_permutation.
 */
typedef struct {
  double *x;
  size_t *P; /**< @brief the permutation matrices*/
} transform_vars_permutation_t;

void transform_vars_permutation_evaluate(coco_problem_t *problem, const double *x, double *y) {
  size_t i;
  transform_vars_permutation_t *data;
  coco_problem_t *inner_problem;

  data = (transform_vars_permutation_t *)coco_problem_transformed_get_data(problem);
  inner_problem = coco_problem_transformed_get_inner_problem(problem);
  for (i = 0; i < inner_problem->number_of_variables; ++i) {
    data->x[i] = x[data->P[i]];
  }

  coco_evaluate_function(inner_problem, data->x, y);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

void transform_vars_permutation_free(void *thing) {
  transform_vars_permutation_t *data = (transform_vars_permutation_t *)thing;
  coco_free_memory(data->x);
  coco_free_memory(data->P);
}

coco_problem_t *transform_vars_permutation(coco_problem_t *inner_problem, const size_t *P,
                                                  const size_t number_of_variables) {
  coco_problem_t *problem;
  transform_vars_permutation_t *data;

  assert(number_of_variables == inner_problem->number_of_variables); /*tmp*/

  data = (transform_vars_permutation_t *)coco_allocate_memory(sizeof(*data));
  data->x = coco_allocate_vector(number_of_variables);
  data->P = coco_duplicate_size_t_vector(P, number_of_variables);

  problem = coco_problem_transformed_allocate(inner_problem, data, transform_vars_permutation_free,
                                              "transform_vars_permutation");
  problem->evaluate_function = transform_vars_permutation_evaluate;
  return problem;
}

void transform_vars_inverse_permutation_evaluate(coco_problem_t *problem, const double *x, double *y) {
  size_t i;
  transform_vars_permutation_t *data;
  coco_problem_t *inner_problem;

  data = (transform_vars_permutation_t *)coco_problem_transformed_get_data(problem);
  inner_problem = coco_problem_transformed_get_inner_problem(problem);
  for (i = 0; i < inner_problem->number_of_variables; ++i) {
    data->x[data->P[i]] = x[i];
  }

  coco_evaluate_function(inner_problem, data->x, y);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

coco_problem_t *transform_vars_inverse_permutation(coco_problem_t *inner_problem, const size_t *P,
                                                          const size_t number_of_variables) {
  coco_problem_t *problem;
  transform_vars_permutation_t *data;

  assert(number_of_variables == inner_problem->number_of_variables);

  data = (transform_vars_permutation_t *)coco_allocate_memory(sizeof(*data));
  data->x = coco_allocate_vector(number_of_variables);
  data->P = coco_duplicate_size_t_vector(P, number_of_variables);

  problem = coco_problem_transformed_allocate(inner_problem, data, transform_vars_permutation_free,
                                              "transform_vars_inverse_permutation");
  problem->evaluate_function = transform_vars_inverse_permutation_evaluate;
  return problem;
}
