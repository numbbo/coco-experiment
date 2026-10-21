/**
 * @file f_katsuura.cpp
 * @brief Implementation of the Katsuura function and problem.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "coco_utilities.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_norm_by_dim.cpp"
#include "transform_obj_penalize.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_affine.cpp"
#include "transform_vars_shift.cpp"

/**
 * @brief Implements the Katsuura function without connections to any COCO
 * structures.
 */
static double f_katsuura_raw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  const double exponent = 10. / pow(static_cast<double>(number_of_variables), 1.2);

  /* Computation core */
  double result = 1.0;
  for (size_t i = 0; i < number_of_variables; ++i) {
    double tmp = 0;
    for (int j = 1; j < 33; ++j) {
      const double tmp2 = ldexp(1.0, j); /* 2^j */
      tmp += fabs(tmp2 * x[i] - coco_double_round(tmp2 * x[i])) / tmp2;
    }
    tmp = 1.0 + (static_cast<double>(i) + 1) * tmp;
    result *= pow(tmp, exponent);
  }
  result =
      10. / (static_cast<double>(number_of_variables)) / (static_cast<double>(number_of_variables)) * (-1. + result);

  return result;
}

/**
 * @brief The Katsuura function.
 */
class f_katsuura_t final : public coco_problem_s {
public:
  explicit f_katsuura_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("Katsuura function", -5.0, 5.0, 0.0);
    set_id("%s_d%02zu", "katsuura", number_of_variables);

    /* Compute best solution */
    evaluate_function(best_parameter, best_value);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives == 1);
    y[0] = f_katsuura_raw(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Allocates the basic Katsuura problem.
 */
static coco_problem_ptr_t f_katsuura_allocate(const size_t number_of_variables) {
  return std::make_unique<f_katsuura_t>(number_of_variables);
}

/**
 * @brief Creates the BBOB Katsuura problem.
 */
static coco_problem_ptr_t f_katsuura_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                           const size_t instance, const long rseed,
                                                           const coco_bbob_variant_t variant,
                                                           const char *problem_id_template,
                                                           const char *problem_name_template) {
  coco_problem_ptr_t problem;
  coco_matrix_t M(dimension, dimension);
  std::vector<double> b(dimension);

  const double penalty_factor = 1.0;

  std::vector<double> xopt(dimension);
  double fopt = bbob2009_compute_fopt(function, instance);
  bbob_compute_xopt(xopt, rseed, variant);

  const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
  const coco_matrix_t rot2 = bbob2009_compute_rotation(rseed, dimension);

  std::vector<double> scales(dimension);
  for (size_t k = 0; k < dimension; ++k) {
    const double exponent = 1.0 * static_cast<int>(k) / (static_cast<double>(dimension) - 1.0);
    scales[k] = pow(sqrt(100), exponent);
  }
  M = coco_matrix_product_with_diagonal(rot1, scales, rot2);

  problem = f_katsuura_allocate(dimension);
  problem = transform_obj_shift(std::move(problem), fopt); /*There is no shift 'fopt' in the definition*/
  problem = transform_vars_affine(std::move(problem), M, b);
  problem = transform_vars_shift(std::move(problem), xopt, 0);
  problem = transform_obj_penalize(std::move(problem), penalty_factor);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("5-weakly-structured");

  return problem;
}

/**
 * @brief Creates the BBOB permuted block-rotated Katsuura problem.
 */
static coco_problem_ptr_t f_katsuura_permblockdiag_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                         const size_t instance, const long rseed,
                                                                         const coco_bbob_variant_t variant,
                                                                         const char *problem_id_template,
                                                                         const char *problem_name_template) {
  coco_problem_ptr_t problem;
  const double penalty_factor = 1.0;
  std::vector<size_t> P11(dimension);
  std::vector<size_t> P12(dimension);
  std::vector<size_t> P21(dimension);
  std::vector<size_t> P22(dimension);
  size_t nb_blocks1;
  size_t nb_blocks2;

  std::vector<size_t> block_sizes1 = coco_get_block_sizes(&nb_blocks1, dimension, "bbob-largescale");
  std::vector<size_t> block_sizes2 = coco_get_block_sizes(&nb_blocks2, dimension, "bbob-largescale");
  size_t swap_range1 = coco_get_swap_range(dimension, "bbob-largescale");
  size_t swap_range2 = coco_get_swap_range(dimension, "bbob-largescale");
  size_t nb_swaps1 = coco_get_nb_swaps(dimension, "bbob-largescale");
  size_t nb_swaps2 = coco_get_nb_swaps(dimension, "bbob-largescale");

  std::vector<double> xopt(dimension);
  double fopt = bbob2009_compute_fopt(function, instance);
  bbob_compute_xopt(xopt, rseed, variant);

  const coco_block_diagonal_matrix_t B1 = coco_compute_blockrotation(rseed + 1000000, block_sizes1);
  const coco_block_diagonal_matrix_t B2 = coco_compute_blockrotation(rseed, block_sizes2);

  coco_compute_truncated_uniform_swap_permutation(P11, rseed + 3000000, nb_swaps1, swap_range1);
  coco_compute_truncated_uniform_swap_permutation(P12, rseed + 4000000, nb_swaps1, swap_range1);
  coco_compute_truncated_uniform_swap_permutation(P21, rseed + 5000000, nb_swaps2, swap_range2);
  coco_compute_truncated_uniform_swap_permutation(P22, rseed + 6000000, nb_swaps2, swap_range2);

  problem = f_katsuura_allocate(dimension);
  problem = transform_vars_permutation(std::move(problem), P22);
  problem = transform_vars_blockrotation(std::move(problem), B1);
  problem = transform_vars_permutation(std::move(problem), P21);
  problem = transform_vars_conditioning(std::move(problem), 100.0);
  problem = transform_vars_permutation(std::move(problem), P12);
  problem = transform_vars_blockrotation(std::move(problem), B2);
  problem = transform_vars_permutation(std::move(problem), P11);
  problem = transform_vars_shift(std::move(problem), xopt, 0);

  problem = transform_obj_penalize(std::move(problem), penalty_factor);
  problem = transform_obj_shift(std::move(problem),
                                fopt); /*TODO: documentation, there is no fopt in the definition of this function*/

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("5-weakly-structured");

  return problem;
}
