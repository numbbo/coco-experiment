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
 * @brief The generalized discus function.
 */
class f_discus_generalized_t final : public coco_problem_s {
public:
  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double raw(const double *x) const;

  f_discus_generalized_t(const size_t number_of_variables, size_t proportion_short_axes_denom)
      : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("generalized discus function", -5.0, 5.0, 0.0);
    set_id("%s_d%04zu", "discus_generalized", number_of_variables);
    this->proportion_short_axes_denom = proportion_short_axes_denom;

    /* Compute best solution */
    evaluate_function(best_parameter.data(), best_value.data());
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = raw(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

private:
  size_t proportion_short_axes_denom{};
};

/**
 * @brief Implements the generalized discus function without connections to any
 * COCO structures.
 */
double f_discus_generalized_t::raw(const double *x) const {

  static const double condition = 1.0e6;
  double result = 0;
  size_t nb_short_axes = number_of_variables / proportion_short_axes_denom;
  if (number_of_variables % proportion_short_axes_denom != 0) {
    nb_short_axes += 1;
  }

  for (size_t i = 0; i < nb_short_axes; ++i) {
    result += x[i] * x[i];
  }
  result *= condition;

  for (size_t i = nb_short_axes; i < number_of_variables; ++i) {
    result += x[i] * x[i];
  }

  return result;
}

/**
 * @brief Allocates the basic generalized discus problem.
 */
static coco_problem_ptr_t f_discus_generalized_allocate(const size_t number_of_variables,
                                                        size_t proportion_short_axes_denom) {
  return std::make_unique<f_discus_generalized_t>(number_of_variables, proportion_short_axes_denom);
}

/**
 * @brief Creates the BBOB generalized permuted block-rotated discus problem.
 */
static coco_problem_ptr_t f_discus_generalized_permblockdiag_bbob_problem_allocate(
    const size_t function, const size_t dimension, const size_t instance, const long rseed,
    const coco_bbob_variant_t variant, const char *problem_id_template, const char *problem_name_template) {
  coco_problem_ptr_t problem;
  std::vector<size_t> P1(dimension);
  std::vector<size_t> P2(dimension);
  size_t nb_blocks;

  const size_t proportion_short_axes_denom = 40;

  std::vector<size_t> block_sizes = coco_get_block_sizes(&nb_blocks, dimension, "bbob-largescale");
  size_t swap_range = coco_get_swap_range(dimension, "bbob-largescale");
  size_t nb_swaps = coco_get_nb_swaps(dimension, "bbob-largescale");

  std::vector<double> xopt(dimension);
  bbob_compute_xopt(xopt.data(), rseed, dimension, variant);
  double fopt = bbob2009_compute_fopt(function, instance);

  const coco_block_diagonal_matrix_t B = coco_compute_blockrotation(rseed + 1000000, block_sizes.data(), nb_blocks);
  coco_compute_truncated_uniform_swap_permutation(
      P1.data(), rseed + 2000000, dimension, nb_swaps, swap_range);
  coco_compute_truncated_uniform_swap_permutation(
      P2.data(), rseed + 3000000, dimension, nb_swaps, swap_range);

  problem =
      f_discus_generalized_allocate(dimension, proportion_short_axes_denom);
  problem = transform_vars_oscillate(std::move(problem));
  problem = transform_vars_permutation(std::move(problem), P2.data(), dimension);
  problem = transform_vars_blockrotation(std::move(problem), B);
  problem = transform_vars_permutation(std::move(problem), P1.data(), dimension);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);

  problem = transform_obj_norm_by_dim(std::move(problem));
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("3-ill-conditioned");

  return problem;
}
