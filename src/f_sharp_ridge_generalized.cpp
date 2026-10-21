/**
 * @file f_sharp_ridge_generalized.cpp
 * @brief Implementation of the generalized sharp ridge function and problem.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_norm_by_dim.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_affine.cpp"
#include "transform_vars_blockrotation.cpp"
#include "transform_vars_conditioning.cpp"
#include "transform_vars_permutation.cpp"
#include "transform_vars_shift.cpp"

/**
 * @brief The sharp ridge function.
 */
class f_sharp_ridge_generalized_t final : public coco_problem_s {
public:
  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double raw(std::span<const double> x) const;

  f_sharp_ridge_generalized_t(const size_t number_of_variables, size_t proportion_of_linear_dims)
      : coco_problem_s(number_of_variables, 1, 0) {
    /* Wassim: proportion_of_linear_dims should probably be allowed to be non-integer */
    init_from_scalars("sharp ridge function", -5.0, 5.0, 0.0);
    set_id("%s_d%02zu", "sharp_ridge_generalized", number_of_variables);
    this->proportion_of_linear_dims = proportion_of_linear_dims;

    /* Compute best solution */
    evaluate_function(best_parameter, best_value);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives == 1);
    y[0] = raw(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

private:
  size_t proportion_of_linear_dims{};
};

/**
 * @brief Implements the generalized sharp ridge function without connections to
 * any COCO structures.
 */
double f_sharp_ridge_generalized_t::raw(std::span<const double> x) const {

  static const double alpha = 100.0;
  size_t i = 0, number_linear_dimensions;

  assert(number_of_variables > 1);

  double result = 0.0;
  number_linear_dimensions = number_of_variables / proportion_of_linear_dims;
  if (number_of_variables % proportion_of_linear_dims != 0) {
    number_linear_dimensions += 1;
  }
  for (i = number_linear_dimensions; i < number_of_variables; ++i) {
    result += x[i] * x[i];
  }
  result = alpha * sqrt(result);
  for (i = 0; i < number_linear_dimensions; i++) {
    result += x[i] * x[i];
  }

  return result;
}

/**
 * @brief Allocates the basic sharp ridge problem.
 */
static coco_problem_ptr_t f_sharp_ridge_generalized_allocate(const size_t number_of_variables,
                                                             size_t proportion_of_linear_dims) {
  return std::make_unique<f_sharp_ridge_generalized_t>(number_of_variables, proportion_of_linear_dims);
}

/**
 * @brief Creates the BBOB permuted block-rotated generalized sharp ridge
 * problem
 */
static coco_problem_ptr_t f_sharp_ridge_generalized_permblockdiag_bbob_problem_allocate(
    const size_t function, const size_t dimension, const size_t instance, const long rseed,
    const coco_bbob_variant_t variant, const char *problem_id_template, const char *problem_name_template) {
  coco_problem_ptr_t problem;
  std::vector<size_t> P11(dimension);
  std::vector<size_t> P21(dimension);
  std::vector<size_t> P12(dimension);
  std::vector<size_t> P22(dimension);
  size_t nb_blocks1;
  size_t nb_blocks2;

  const size_t proportion_of_linear_dims = 40;

  std::vector<size_t> block_sizes1 = coco_get_block_sizes(&nb_blocks1, dimension, "bbob-largescale");
  std::vector<size_t> block_sizes2 = coco_get_block_sizes(&nb_blocks2, dimension, "bbob-largescale");
  size_t swap_range = coco_get_swap_range(dimension, "bbob-largescale");
  size_t nb_swaps = coco_get_nb_swaps(dimension, "bbob-largescale");

  std::vector<double> xopt(dimension);
  double fopt = bbob2009_compute_fopt(function, instance);
  bbob_compute_xopt(xopt, rseed, variant);

  const coco_block_diagonal_matrix_t B1 = coco_compute_blockrotation(rseed + 1000000, block_sizes1);
  const coco_block_diagonal_matrix_t B2 = coco_compute_blockrotation(rseed, block_sizes2);

  coco_compute_truncated_uniform_swap_permutation(P11, rseed + 3000000, nb_swaps, swap_range);
  coco_compute_truncated_uniform_swap_permutation(P21, rseed + 4000000, nb_swaps, swap_range);
  coco_compute_truncated_uniform_swap_permutation(P12, rseed + 5000000, nb_swaps, swap_range);
  coco_compute_truncated_uniform_swap_permutation(P22, rseed + 6000000, nb_swaps, swap_range);

  problem = f_sharp_ridge_generalized_allocate(dimension, proportion_of_linear_dims);
  problem = transform_vars_permutation(std::move(problem), P21); /* LIFO */
  problem = transform_vars_blockrotation(std::move(problem), B1);
  problem = transform_vars_permutation(std::move(problem), P11);
  problem = transform_vars_conditioning(std::move(problem), 10.0);
  problem =
      transform_vars_permutation(std::move(problem), P22); /*Consider replacing P11 and 22 by a single permutation P3*/
  problem = transform_vars_blockrotation(std::move(problem), B2);
  problem = transform_vars_permutation(std::move(problem), P12);
  problem = transform_vars_shift(std::move(problem), xopt, 0);

  problem = transform_obj_norm_by_dim(std::move(problem));
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("3-ill-conditioned");

  return problem;
}
