/**
 * @file f_schaffers.cpp
 * @brief Implementation of the Schaffer's F7 function and problem,
 * transformations not implemented for the moment.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_norm_by_dim.cpp"
#include "transform_obj_penalize.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_affine.cpp"
#include "transform_vars_asymmetric.cpp"
#include "transform_vars_blockrotation.cpp"

/**
 * @brief Implements the Schaffer's F7 function without connections to any COCO
 * structures.
 */
static double f_schaffers_raw(const double *x,
                              const size_t number_of_variables) {

  size_t i = 0;

  assert(number_of_variables > 1);

  /* Computation core */
  double result = 0.0;
  for (i = 0; i < number_of_variables - 1; ++i) {
    const double tmp = x[i] * x[i] + x[i + 1] * x[i + 1];
    if (coco_is_inf(tmp) && coco_is_nan(sin(50.0 * pow(tmp, 0.1)))) /* sin(inf) -> nan */
      /* the second condition is necessary to pass the integration tests under Windows and Linux */
      return tmp;
    result += pow(tmp, 0.25) * (1.0 + pow(sin(50.0 * pow(tmp, 0.1)), 2.0));
  }
  result = pow(result / (static_cast<double>(number_of_variables) - 1.0), 2.0);

  return result;
}

/**
 * @brief The Schaffer's function.
 */
class f_schaffers_t final : public coco_problem_s {
public:
  explicit f_schaffers_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("Schaffer's function", -5.0, 5.0, 0.0);
    set_id("%s_d%02zu", "schaffers", number_of_variables);

    /* Compute best solution */
    evaluate_function(best_parameter.data(), best_value.data());
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_schaffers_raw(x, number_of_variables);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Allocates the basic Schaffer's F7 problem.
 */
static coco_problem_ptr_t f_schaffers_allocate(const size_t number_of_variables) {
  return std::make_unique<f_schaffers_t>(number_of_variables);
}

/**
 * @brief Creates the BBOB Schaffer's F7 problem.
 */
static coco_problem_ptr_t f_schaffers_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                            const size_t instance, const long rseed,
                                                            const coco_bbob_variant_t variant,
                                                            const double conditioning, const double penalty_scale,
                                                            const char *problem_id_template,
                                                            const char *problem_name_template) {
  coco_problem_ptr_t problem;
  coco_matrix_t M(dimension, dimension);
  std::vector<double> b(dimension);


  std::vector<double> xopt(dimension);
  double fopt = bbob2009_compute_fopt(function, instance);
  bbob_compute_xopt(xopt.data(), rseed, dimension, variant);

  const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
  const coco_matrix_t rot2 = bbob2009_compute_rotation(rseed, dimension);
  for (size_t i = 0; i < dimension; ++i) {
    b[i] = 0.0;
    for (size_t j = 0; j < dimension; ++j) {
      double exponent = 1.0 * static_cast<int>(i) / (static_cast<double>(dimension) - 1.0);
      M(i, j) = rot2(i, j) * pow(sqrt(conditioning), exponent);
    }
  }

  problem = f_schaffers_allocate(dimension);
  problem = transform_obj_shift(std::move(problem), fopt);
  problem = transform_vars_affine(std::move(problem), M, b);
  problem = transform_vars_asymmetric(std::move(problem), 0.5);
  M = rot1;
  problem = transform_vars_affine(std::move(problem), M, b);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);
  problem = transform_obj_penalize(std::move(problem), penalty_scale);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("4-multi-modal");

  return problem;
}

/**
 * @brief Creates the BBOB permuted block-rotated Schaffer's F7 problem.
 */
static coco_problem_ptr_t f_schaffers_permblockdiag_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                          const size_t instance, const long rseed,
                                                                          const coco_bbob_variant_t variant,
                                                                          const double conditioning,
                                                                          const char *problem_id_template,
                                                                          const char *problem_name_template) {
  coco_problem_ptr_t problem;
  std::vector<size_t> P11(dimension);
  std::vector<size_t> P21(dimension);
  std::vector<size_t> P12(dimension);
  std::vector<size_t> P22(dimension);
  size_t nb_blocks1;
  size_t nb_blocks2;

  const double penalty_factor = 10.0;

  std::vector<size_t> block_sizes1 = coco_get_block_sizes(&nb_blocks1, dimension, "bbob-largescale");
  std::vector<size_t> block_sizes2 = coco_get_block_sizes(&nb_blocks2, dimension, "bbob-largescale");
  size_t swap_range = coco_get_swap_range(dimension, "bbob-largescale");
  size_t nb_swaps = coco_get_nb_swaps(dimension, "bbob-largescale");

  std::vector<double> xopt(dimension);
  double fopt = bbob2009_compute_fopt(function, instance);
  bbob_compute_xopt(xopt.data(), rseed, dimension, variant);

  const coco_block_diagonal_matrix_t B1 = coco_compute_blockrotation(rseed + 1000000, block_sizes1.data(), nb_blocks1);
  const coco_block_diagonal_matrix_t B2 = coco_compute_blockrotation(rseed, block_sizes2.data(), nb_blocks2);

  coco_compute_truncated_uniform_swap_permutation(P11.data(), rseed + 3000000, dimension, nb_swaps, swap_range);
  coco_compute_truncated_uniform_swap_permutation(P21.data(), rseed + 4000000, dimension, nb_swaps, swap_range);
  coco_compute_truncated_uniform_swap_permutation(P12.data(), rseed + 5000000, dimension, nb_swaps, swap_range);
  coco_compute_truncated_uniform_swap_permutation(P22.data(), rseed + 6000000, dimension, nb_swaps, swap_range);

  problem = f_schaffers_allocate(dimension);
  problem = transform_vars_conditioning(std::move(problem), conditioning);
  problem = transform_vars_permutation(std::move(problem), P21.data(), dimension);
  problem = transform_vars_blockrotation(std::move(problem), B2);
  problem = transform_vars_permutation(std::move(problem), P11.data(), dimension);

  problem = transform_vars_asymmetric(std::move(problem), 0.5);
  problem = transform_vars_permutation(std::move(problem), P22.data(), dimension);
  problem = transform_vars_blockrotation(std::move(problem), B1);
  problem = transform_vars_permutation(std::move(problem), P12.data(), dimension);

  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);
  problem = transform_obj_penalize(std::move(problem), penalty_factor);
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("4-multi-modal");

  return problem;
}
