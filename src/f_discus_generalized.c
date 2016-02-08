/**
 * @file f_discus.c
 * @brief Implementation of the discus function and problem.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.c"
#include "suite_bbob_legacy_code.c"
#include "transform_vars_oscillate.c"
#include "transform_vars_shift.c"
#include "transform_obj_shift.c"

#include "transform_vars_permblockdiag.c"
#include "transform_obj_norm_by_dim.c"

#define proportion_short_axes_denom 40 /* Make it a parameter of f_discus_generalized_raw? */

/**
 * @brief Implements the generalized discus function without connections to any COCO structures.
 */
static double f_discus_generalized_raw(const double *x, const size_t number_of_variables) {

  static const double condition = 1.0e6;
  size_t i, nb_short_axes;
  double result;
  result = 0;
  nb_short_axes = number_of_variables / proportion_short_axes_denom;
  if (number_of_variables % proportion_short_axes_denom != 0) {
    nb_short_axes += 1;
  }

  for (i = 0; i < nb_short_axes; ++i) {
    result += x[i] * x[i];
  }
  result *= condition;

  for (i = nb_short_axes; i < number_of_variables; ++i) {
    result += x[i] * x[i];
  }

  return result;
}

/**
 * @brief Uses the generalized raw function to evaluate the COCO problem.
 */
static void f_discus_generalized_evaluate(coco_problem_t *problem, const double *x, double *y) {
  assert(problem->number_of_objectives == 1);
  y[0] = f_discus_generalized_raw(x, problem->number_of_variables);
  assert(y[0] + 1e-13 >= problem->best_value[0]);
}

/**
 * @brief Allocates the basic generalized discus problem.
 */
static coco_problem_t *f_discus_generalized_allocate(const size_t number_of_variables) {

  coco_problem_t *problem = coco_problem_allocate_from_scalars("generalized discus function",
                                                               f_discus_generalized_evaluate, NULL, number_of_variables, -5.0, 5.0, 0.0);
  coco_problem_set_id(problem, "%s_d%02lu", "discus_generalized", number_of_variables);

  /* Compute best solution */
  f_discus_generalized_evaluate(problem, problem->best_parameter, problem->best_value);
  return problem;
}



/**
 * @brief Creates the BBOB generalized discus problem.
 */
static coco_problem_t *f_discus_generalized_permblockdiag_bbob_problem_allocate(const size_t function,
                                                                                const size_t dimension,
                                                                                const size_t instance,
                                                                                const long rseed,
                                                                                const char *problem_id_template,
                                                                                const char *problem_name_template) {
  double *xopt, fopt;
  coco_problem_t *problem = NULL;
  double **B;
  const double *const *B_copy;
  size_t *P1 = coco_allocate_vector_size_t(dimension);
  size_t *P2 = coco_allocate_vector_size_t(dimension);
  size_t *block_sizes;
  size_t nb_blocks;
  size_t swap_range;
  size_t nb_swaps;

  block_sizes = ls_get_block_sizes(&nb_blocks, dimension);
  swap_range = ls_get_swap_range(dimension);
  nb_swaps = ls_get_nb_swaps(dimension);

  xopt = coco_allocate_vector(dimension);
  bbob2009_compute_xopt(xopt, rseed, dimension);
  fopt = bbob2009_compute_fopt(function, instance);

  B = ls_allocate_blockmatrix(dimension, block_sizes, nb_blocks);
  B_copy = (const double *const *)B;

  ls_compute_blockrotation(B, rseed + 1000000, dimension, block_sizes, nb_blocks);
  ls_compute_truncated_uniform_swap_permutation(P1, rseed + 2000000, dimension, nb_swaps, swap_range);
  ls_compute_truncated_uniform_swap_permutation(P2, rseed + 3000000, dimension, nb_swaps, swap_range);

  problem = f_discus_generalized_allocate(dimension);
  problem = transform_vars_oscillate(problem);
  problem = transform_vars_permblockdiag(problem, B_copy, P1, P2, dimension, block_sizes, nb_blocks);
  problem = transform_vars_shift(problem, xopt, 0);

  problem = transform_obj_norm_by_dim(problem);
  problem = transform_obj_shift(problem, fopt);

  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "large_scale_block_rotated");/*TODO: no large scale prefix*/

  ls_free_block_matrix(B, dimension);
  coco_free_memory(P1);
  coco_free_memory(P2);
  coco_free_memory(block_sizes);
  coco_free_memory(xopt);
  return problem;
}

