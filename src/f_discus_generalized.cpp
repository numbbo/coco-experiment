/**
 * @file f_discus_generalized.cpp
 * @brief Implementation of the discus function and problem.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_norm_by_dim.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_blockrotation.cpp"
#include "transform_vars_oscillate.cpp"
#include "transform_vars_permutation.cpp"
#include "transform_vars_shift.cpp"

/**
 * @brief Parameters of the generalized discus function.
 */
struct f_discus_generalized_data_t {
  size_t proportion_short_axes_denom;
};

/**
 * @brief Implements the generalized discus function without connections to any
 * COCO structures.
 */
static double f_discus_generalized_raw(const double *x, const size_t number_of_variables,
                                       f_discus_generalized_data_t *data) {

  static const double condition = 1.0e6;
  size_t i, nb_short_axes;
  double result;
  result = 0;
  nb_short_axes =
      number_of_variables /
      data->proportion_short_axes_denom;
  if (number_of_variables %
          data->proportion_short_axes_denom !=
      0) {
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
 * @brief The generalized discus function.
 */
struct f_discus_generalized_t final : public coco_problem_s {
  f_discus_generalized_data_t params;

  f_discus_generalized_t(const size_t number_of_variables, size_t proportion_short_axes_denom)
      : coco_problem_s(number_of_variables, 1, 0) {
    coco_problem_init_from_scalars(this, "generalized discus function", -5.0, 5.0, 0.0);
    coco_problem_set_id(this, "%s_d%04lu", "discus_generalized", number_of_variables);
    params.proportion_short_axes_denom = proportion_short_axes_denom;

    /* Compute best solution */
    evaluate_function(best_parameter.data(), best_value.data());
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_discus_generalized_raw(x, number_of_variables, &params);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Allocates the basic generalized discus problem.
 */
static coco_problem_t *f_discus_generalized_allocate(const size_t number_of_variables,
                                                     size_t proportion_short_axes_denom) {
  return new f_discus_generalized_t(number_of_variables, proportion_short_axes_denom);
}

/**
 * @brief Creates the BBOB generalized permuted block-rotated discus problem.
 */
static coco_problem_t *f_discus_generalized_permblockdiag_bbob_problem_allocate(const size_t function,
                                                                                const size_t dimension,
                                                                                const size_t instance,
                                                                                const long rseed,
                                                                                const char *problem_id_template,
                                                                                const char *problem_name_template) {
  double fopt;
  coco_problem_t *problem = nullptr;
  std::vector<size_t> P1(dimension);
  std::vector<size_t> P2(dimension);
  size_t nb_blocks;
  size_t swap_range;
  size_t nb_swaps;

  const size_t proportion_short_axes_denom = 40;

  std::vector<size_t> block_sizes = coco_get_block_sizes(&nb_blocks, dimension, "bbob-largescale");
  swap_range = coco_get_swap_range(dimension, "bbob-largescale");
  nb_swaps = coco_get_nb_swaps(dimension, "bbob-largescale");

  std::vector<double> xopt(dimension);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt.data(), rseed, dimension);
  } else {
    bbob2009_compute_xopt(xopt.data(), rseed, dimension);
  }
  fopt = bbob2009_compute_fopt(function, instance);

  coco_matrix_t B = coco_block_matrix(dimension, block_sizes.data(), nb_blocks);

  coco_compute_blockrotation(B.rows(), rseed + 1000000, dimension, block_sizes.data(),
                             nb_blocks);
  coco_compute_truncated_uniform_swap_permutation(
      P1.data(), rseed + 2000000, dimension, nb_swaps, swap_range);
  coco_compute_truncated_uniform_swap_permutation(
      P2.data(), rseed + 3000000, dimension, nb_swaps, swap_range);

  problem =
      f_discus_generalized_allocate(dimension, proportion_short_axes_denom);
  problem = transform_vars_oscillate(problem);
  problem = transform_vars_permutation(problem, P2.data(), dimension);
  problem = transform_vars_blockrotation(problem, B.rows(), dimension,
                                         block_sizes.data(), nb_blocks);
  problem = transform_vars_permutation(problem, P1.data(), dimension);
  problem = transform_vars_shift(problem, xopt.data(), 0);

  problem = transform_obj_norm_by_dim(problem);
  problem = transform_obj_shift(problem, fopt);

  coco_problem_set_id(problem, problem_id_template, function, instance,
                      dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance,
                        dimension);
  coco_problem_set_type(problem, "3-ill-conditioned");

  return problem;
}
