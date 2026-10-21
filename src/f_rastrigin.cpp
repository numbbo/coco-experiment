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
  double result;
  double sum1 = 0.0, sum2 = 0.0;

  if (coco_vector_contains_nan(x, number_of_variables))
    return NAN;

  for (i = 0; i < number_of_variables; ++i) {
    sum1 += cos(coco_two_pi * x[i]);
    sum2 += x[i] * x[i];
  }
  if (coco_is_inf(sum2)) /* cos(inf) -> nan */
    return sum2;
  result = 10.0 * ((double)(long)number_of_variables - sum1) + sum2;

  return result;
}

/**
 * @brief The Rastrigin function.
 */
struct f_rastrigin_t final : public coco_problem_s {
  explicit f_rastrigin_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    coco_problem_init_from_scalars(this, "Rastrigin function", -5.0, 5.0, 0.0);
    /* TODO: make sure the gradient is computed correctly for the rotated Rastrigin */
    coco_problem_set_id(this, "%s_d%02lu", "rastrigin", number_of_variables);

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

    size_t i;

    for (i = 0; i < number_of_variables; ++i) {
      y[i] = 2.0 * (10. * coco_pi * sin(coco_two_pi * x[i]) + x[i]);
    }
  }
};

/**
 * @brief Allocates the basic Rastrigin problem.
 */
static coco_problem_t *f_rastrigin_allocate(const size_t number_of_variables) {
  return new f_rastrigin_t(number_of_variables);
}

/**
 * @brief Creates the BBOB Rastrigin problem.
 */
static coco_problem_t *f_rastrigin_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                         const size_t instance, const long rseed,
                                                         const char *problem_id_template,
                                                         const char *problem_name_template) {

  double fopt;
  coco_problem_t *problem = nullptr;

  std::vector<double> xopt(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt.data(), rseed, dimension);
  } else {
    bbob2009_compute_xopt(xopt.data(), rseed, dimension);
  }

  problem = f_rastrigin_allocate(dimension);
  problem = transform_vars_conditioning(problem, 10.0);
  problem = transform_vars_asymmetric(problem, 0.2);
  problem = transform_vars_oscillate(problem);
  problem = transform_vars_shift(problem, xopt.data(), 0);

  /*if large scale test-bed, normalize by dim*/
  if (coco_strfind(problem_name_template, "BBOB large-scale suite") >= 0) {
    problem = transform_obj_norm_by_dim(problem);
  }
  problem = transform_obj_shift(problem, fopt);

  coco_problem_set_id(problem, problem_id_template, function, instance,
                      dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance,
                        dimension);
  coco_problem_set_type(problem, "1-separable");

  return problem;
}

/**
 * @brief Creates the BBOB rotated Rastrigin problem.
 */
static coco_problem_t *f_rastrigin_rotated_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                 const size_t instance, const long rseed,
                                                                 const char *problem_id_template,
                                                                 const char *problem_name_template) {
  double fopt;
  coco_problem_t *problem = nullptr;
  size_t i, j, k;
  std::vector<double> M(dimension * dimension);
  std::vector<double> b(dimension);
  double *current_row;

  std::vector<double> xopt(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt.data(), rseed, dimension);
  } else {
    bbob2009_compute_xopt(xopt.data(), rseed, dimension);
  }

  coco_matrix_t rot1(dimension, dimension);
  coco_matrix_t rot2(dimension, dimension);
  bbob2009_compute_rotation(rot1.rows(), rseed + 1000000, dimension);
  bbob2009_compute_rotation(rot2.rows(), rseed, dimension);
  for (i = 0; i < dimension; ++i) {
    b[i] = 0.0;
    current_row = M.data() + i * dimension;
    for (j = 0; j < dimension; ++j) {
      current_row[j] = 0.0;
      for (k = 0; k < dimension; ++k) {
        double exponent = 1.0 * (int)k / ((double)(long)dimension - 1.0);
        current_row[j] += rot1[i][k] * pow(sqrt(10), exponent) * rot2[k][j];
      }
    }
  }

  problem = f_rastrigin_allocate(dimension);
  problem = transform_obj_shift(problem, fopt);
  problem = transform_vars_affine(problem, M.data(), b.data(), dimension);
  problem = transform_vars_asymmetric(problem, 0.2);
  problem = transform_vars_oscillate(problem);
  bbob2009_copy_rotation_matrix(rot1.rows(), M.data(), b.data(), dimension);
  problem = transform_vars_affine(problem, M.data(), b.data(), dimension);
  problem = transform_vars_shift(problem, xopt.data(), 0);


  coco_problem_set_id(problem, problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem, "4-multi-modal");

  return problem;
}

/**
 * @brief Creates the BBOB rotated Rastrigin problem for large scale.
 */
static coco_problem_t *f_rastrigin_permblockdiag_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                       const size_t instance, const long rseed,
                                                                       const char *problem_id_template,
                                                                       const char *problem_name_template) {
  double fopt;
  coco_problem_t *problem = nullptr;
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
  fopt = bbob2009_compute_fopt(function, instance);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt.data(), rseed, dimension);
  } else {
    bbob2009_compute_xopt(xopt.data(), rseed, dimension);
  }

  coco_matrix_t B1 = coco_block_matrix(dimension, block_sizes1.data(), nb_blocks1);
  coco_matrix_t B2 = coco_block_matrix(dimension, block_sizes2.data(), nb_blocks2);

  coco_compute_blockrotation(B1.rows(), rseed + 1000000, dimension, block_sizes1.data(),
                             nb_blocks1);
  coco_compute_blockrotation(B2.rows(), rseed, dimension, block_sizes2.data(), nb_blocks2);

  coco_compute_truncated_uniform_swap_permutation(P11.data(), rseed + 3000000, dimension, nb_swaps1, swap_range1);
  coco_compute_truncated_uniform_swap_permutation(P12.data(), rseed + 4000000, dimension, nb_swaps1, swap_range1);
  coco_compute_truncated_uniform_swap_permutation(P21.data(), rseed + 5000000, dimension, nb_swaps2, swap_range2);
  coco_compute_truncated_uniform_swap_permutation(P22.data(), rseed + 6000000, dimension, nb_swaps2, swap_range2);

  problem = f_rastrigin_allocate(dimension);
  problem = transform_vars_permutation(problem, P12.data(), dimension);
  problem = transform_vars_blockrotation(problem, B1.rows(), dimension,
                                         block_sizes1.data(), nb_blocks1);
  problem = transform_vars_permutation(problem, P11.data(), dimension);
  problem = transform_vars_conditioning(problem, 10.0);
  problem = transform_vars_permutation(problem, P22.data(), dimension);
  problem = transform_vars_blockrotation(problem, B2.rows(), dimension,
                                         block_sizes2.data(), nb_blocks2);
  problem = transform_vars_permutation(problem, P21.data(), dimension);
  problem = transform_vars_asymmetric(problem, 0.2);
  problem = transform_vars_oscillate(problem);
  problem = transform_vars_permutation(problem, P12.data(), dimension);
  problem = transform_vars_blockrotation(problem, B1.rows(), dimension,
                                         block_sizes1.data(), nb_blocks1);
  problem = transform_vars_permutation(problem, P11.data(), dimension);
  problem = transform_vars_shift(problem, xopt.data(), 0);

  problem = transform_obj_norm_by_dim(problem);
  problem = transform_obj_shift(problem, fopt);

  coco_problem_set_id(problem, problem_id_template, function, instance,
                      dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance,
                        dimension);
  coco_problem_set_type(problem, "4-multi-modal");

  return problem;
}

/**
 * @brief Computes xopt for constrained Rastrigin (alternative to
 * bbob2009_compute_xopt()) xopt is a vector of dim uniform random integers
 */
static void f_rastrigin_cons_compute_xopt(double *xopt, const long rseed,
                                          const size_t dim) {

  size_t i;

  bbob2009_unif(xopt, dim, rseed);

  for (i = 0; i < dim; ++i) {
    xopt[i] = 10 * xopt[i] - 5;
    xopt[i] = (int)xopt[i];
  }

  /* In case (0, ..., 0) is sampled, set xopt to a different value */
  if (coco_vector_is_zero(xopt, dim))
    for (i = 0; i < dim; ++i) {
      xopt[i] = (int)(i % 9) - 4;
    }
}

/**
 * @brief Creates the Rastrigin problem for the constrained BBOB suite.
 */
static coco_problem_t *f_rastrigin_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                              const size_t instance, const long rseed,
                                                              const char *problem_id_template,
                                                              const char *problem_name_template) {

  double fopt;
  coco_problem_t *problem = nullptr;

  std::vector<double> xopt(dimension);
  f_rastrigin_cons_compute_xopt(xopt.data(), rseed, dimension);
  fopt = bbob2009_compute_fopt(function, instance);

  problem = f_rastrigin_allocate(dimension);
  problem = transform_vars_shift(problem, xopt.data(), 0);
  problem = transform_obj_shift(problem, fopt);

  coco_problem_set_id(problem, problem_id_template, function, instance,
                      dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance,
                        dimension);
  coco_problem_set_type(problem, "1-separable");

  return problem;
}

/**
 * @brief Creates the rotated Rastrigin problem for the constrained BBOB suite.
 */
static coco_problem_t *f_rastrigin_rotated_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                      const size_t instance, const long rseed,
                                                                      const char *problem_id_template,
                                                                      const char *problem_name_template) {

  double fopt;
  coco_problem_t *problem = nullptr;

  std::vector<double> M(dimension * dimension);
  std::vector<double> b(dimension);

  coco_matrix_t rot1(dimension, dimension);
  bbob2009_compute_rotation(rot1.rows(), rseed + 1000000, dimension);
  bbob2009_copy_rotation_matrix(rot1.rows(), M.data(), b.data(), dimension);

  std::vector<double> xopt(dimension);
  f_rastrigin_cons_compute_xopt(xopt.data(), rseed, dimension);
  fopt = bbob2009_compute_fopt(function, instance);

  problem = f_rastrigin_allocate(dimension);
  problem = transform_vars_shift(problem, xopt.data(), 0);
  problem = transform_vars_affine(problem, M.data(), b.data(), dimension); /* Rotate after shifting so R dot xopt is also rotated */
  problem = transform_obj_shift(problem, fopt);

  coco_problem_set_id(problem, problem_id_template, function, instance,
                      dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance,
                        dimension);
  coco_problem_set_type(problem, "1-separable");

  return problem;
}
