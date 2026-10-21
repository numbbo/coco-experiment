/**
 * @file f_rastrigin.cpp
 * @brief Implementation of the Rastrigin function and problem.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "coco_utilities.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_norm_by_dim.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_affine.cpp"
#include "transform_vars_asymmetric.cpp"
#include "transform_vars_blockrotation.cpp"
#include "transform_vars_conditioning.cpp"
#include "transform_vars_oscillate.cpp"
#include "transform_vars_permutation.cpp"
#include "transform_vars_shift.cpp"

/**
 * @brief Implements the Rastrigin function without connections to any COCO
 * structures.
 */
static double f_rastrigin_raw(const double *x,
                              const size_t number_of_variables) {

  size_t i = 0;
  double sum1 = 0.0, sum2 = 0.0;

  for (i = 0; i < number_of_variables; ++i) {
    sum1 += cos(coco_two_pi * x[i]);
    sum2 += x[i] * x[i];
  }
  if (coco_is_inf(sum2)) /* cos(inf) -> nan */
    return sum2;
  double result = 10.0 * (static_cast<double>(number_of_variables) - sum1) + sum2;

  return result;
}

/**
 * @brief The Rastrigin function.
 */
class f_rastrigin_t final : public coco_problem_s {
public:
  explicit f_rastrigin_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("Rastrigin function", -5.0, 5.0, 0.0);
    /* TODO: make sure the gradient is computed correctly for the rotated Rastrigin */
    set_id("%s_d%02zu", "rastrigin", number_of_variables);

    /* Compute best solution */
    evaluate_function(best_parameter.data(), best_value.data());
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_rastrigin_raw(x, number_of_variables);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the gradient of the raw Rastrigin function.
   */
  void evaluate_gradient(const double *x, double *y) override {

    for (size_t i = 0; i < number_of_variables; ++i) {
      y[i] = 2.0 * (10. * coco_pi * sin(coco_two_pi * x[i]) + x[i]);
    }
  }
};

/**
 * @brief Allocates the basic Rastrigin problem.
 */
static coco_problem_ptr_t f_rastrigin_allocate(const size_t number_of_variables) {
  return std::make_unique<f_rastrigin_t>(number_of_variables);
}

/**
 * @brief Creates the BBOB Rastrigin problem.
 */
static coco_problem_ptr_t f_rastrigin_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                            const size_t instance, const long rseed,
                                                            const coco_bbob_variant_t variant,
                                                            const char *problem_id_template,
                                                            const char *problem_name_template) {

  coco_problem_ptr_t problem;

  std::vector<double> xopt(dimension);
  double fopt = bbob2009_compute_fopt(function, instance);
  bbob_compute_xopt(xopt.data(), rseed, dimension, variant);

  problem = f_rastrigin_allocate(dimension);
  problem = transform_vars_conditioning(std::move(problem), 10.0);
  problem = transform_vars_asymmetric(std::move(problem), 0.2);
  problem = transform_vars_oscillate(std::move(problem));
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);

  /*if large scale test-bed, normalize by dim*/
  if (variant == coco_bbob_variant_t::largescale) {
    problem = transform_obj_norm_by_dim(std::move(problem));
  }
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("1-separable");

  return problem;
}

/**
 * @brief Creates the BBOB rotated Rastrigin problem.
 */
static coco_problem_ptr_t f_rastrigin_rotated_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                    const size_t instance, const long rseed,
                                                                    const coco_bbob_variant_t variant,
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
  std::vector<double> scales(dimension);
  for (size_t k = 0; k < dimension; ++k) {
    const double exponent = 1.0 * static_cast<int>(k) / (static_cast<double>(dimension) - 1.0);
    scales[k] = pow(sqrt(10), exponent);
  }
  M = coco_matrix_product_with_diagonal(rot1, scales, rot2);

  problem = f_rastrigin_allocate(dimension);
  problem = transform_obj_shift(std::move(problem), fopt);
  problem = transform_vars_affine(std::move(problem), M, b);
  problem = transform_vars_asymmetric(std::move(problem), 0.2);
  problem = transform_vars_oscillate(std::move(problem));
  M = rot1;
  problem = transform_vars_affine(std::move(problem), M, b);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("4-multi-modal");

  return problem;
}

/**
 * @brief Creates the BBOB rotated Rastrigin problem for large scale.
 */
static coco_problem_ptr_t f_rastrigin_permblockdiag_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                          const size_t instance, const long rseed,
                                                                          const coco_bbob_variant_t variant,
                                                                          const char *problem_id_template,
                                                                          const char *problem_name_template) {
  coco_problem_ptr_t problem;
  std::vector<size_t> P11(dimension);
  std::vector<size_t> P12(dimension);
  std::vector<size_t> P21(dimension);
  std::vector<size_t> P22(dimension);
  size_t nb_blocks1, nb_blocks2;
  size_t swap_range1, swap_range2;
  size_t nb_swaps1, nb_swaps2;

  std::vector<size_t> block_sizes1 = coco_get_block_sizes(&nb_blocks1, dimension, "bbob-largescale");
  std::vector<size_t> block_sizes2 = coco_get_block_sizes(&nb_blocks2, dimension, "bbob-largescale");
  swap_range1 = coco_get_swap_range(dimension, "bbob-largescale");
  swap_range2 = coco_get_swap_range(dimension, "bbob-largescale");
  nb_swaps1 = coco_get_nb_swaps(dimension, "bbob-largescale");
  nb_swaps2 = coco_get_nb_swaps(dimension, "bbob-largescale");

  std::vector<double> xopt(dimension);
  double fopt = bbob2009_compute_fopt(function, instance);
  bbob_compute_xopt(xopt.data(), rseed, dimension, variant);

  const coco_block_diagonal_matrix_t B1 = coco_compute_blockrotation(rseed + 1000000, block_sizes1.data(), nb_blocks1);
  const coco_block_diagonal_matrix_t B2 = coco_compute_blockrotation(rseed, block_sizes2.data(), nb_blocks2);

  coco_compute_truncated_uniform_swap_permutation(P11.data(), rseed + 3000000, dimension, nb_swaps1, swap_range1);
  coco_compute_truncated_uniform_swap_permutation(P12.data(), rseed + 4000000, dimension, nb_swaps1, swap_range1);
  coco_compute_truncated_uniform_swap_permutation(P21.data(), rseed + 5000000, dimension, nb_swaps2, swap_range2);
  coco_compute_truncated_uniform_swap_permutation(P22.data(), rseed + 6000000, dimension, nb_swaps2, swap_range2);

  problem = f_rastrigin_allocate(dimension);
  problem = transform_vars_permutation(std::move(problem), P12.data(), dimension);
  problem = transform_vars_blockrotation(std::move(problem), B1);
  problem = transform_vars_permutation(std::move(problem), P11.data(), dimension);
  problem = transform_vars_conditioning(std::move(problem), 10.0);
  problem = transform_vars_permutation(std::move(problem), P22.data(), dimension);
  problem = transform_vars_blockrotation(std::move(problem), B2);
  problem = transform_vars_permutation(std::move(problem), P21.data(), dimension);
  problem = transform_vars_asymmetric(std::move(problem), 0.2);
  problem = transform_vars_oscillate(std::move(problem));
  problem = transform_vars_permutation(std::move(problem), P12.data(), dimension);
  problem = transform_vars_blockrotation(std::move(problem), B1);
  problem = transform_vars_permutation(std::move(problem), P11.data(), dimension);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);

  problem = transform_obj_norm_by_dim(std::move(problem));
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("4-multi-modal");

  return problem;
}

/**
 * @brief Computes xopt for constrained Rastrigin (alternative to
 * bbob2009_compute_xopt()) xopt is a vector of dim uniform random integers
 */
static void f_rastrigin_cons_compute_xopt(double *xopt, const long rseed,
                                          const size_t dim) {


  bbob2009_unif(xopt, dim, rseed);

  for (size_t i = 0; i < dim; ++i) {
    xopt[i] = 10 * xopt[i] - 5;
    xopt[i] = static_cast<int>(xopt[i]);
  }

  /* In case (0, ..., 0) is sampled, set xopt to a different value */
  if (coco_vector_is_zero(xopt, dim))
    for (size_t i = 0; i < dim; ++i) {
      xopt[i] = static_cast<int>(i % 9) - 4;
    }
}

/**
 * @brief Creates the Rastrigin problem for the constrained BBOB suite.
 */
static coco_problem_ptr_t f_rastrigin_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                 const size_t instance, const long rseed,
                                                                 const char *problem_id_template,
                                                                 const char *problem_name_template) {

  coco_problem_ptr_t problem;

  std::vector<double> xopt(dimension);
  f_rastrigin_cons_compute_xopt(xopt.data(), rseed, dimension);
  double fopt = bbob2009_compute_fopt(function, instance);

  problem = f_rastrigin_allocate(dimension);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("1-separable");

  return problem;
}

/**
 * @brief Creates the rotated Rastrigin problem for the constrained BBOB suite.
 */
static coco_problem_ptr_t f_rastrigin_rotated_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                         const size_t instance, const long rseed,
                                                                         const char *problem_id_template,
                                                                         const char *problem_name_template) {

  coco_problem_ptr_t problem;

  coco_matrix_t M(dimension, dimension);
  std::vector<double> b(dimension);

  const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
  M = rot1;

  std::vector<double> xopt(dimension);
  f_rastrigin_cons_compute_xopt(xopt.data(), rseed, dimension);
  double fopt = bbob2009_compute_fopt(function, instance);

  problem = f_rastrigin_allocate(dimension);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);
  problem = transform_vars_affine(std::move(problem), M, b); /* Rotate after shifting so R dot xopt is also rotated */
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("1-separable");

  return problem;
}
