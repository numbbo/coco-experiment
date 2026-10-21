/**
 * @file f_different_powers.cpp
 * @brief Implementation of the different powers function and problem.
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
#include "transform_vars_permutation.cpp"
#include "transform_vars_shift.cpp"

/**
 * @brief Implements the different powers function without connections to any
 * COCO structures.
 */
static double f_different_powers_raw(const double *x,
                                     const size_t number_of_variables) {

  size_t i;
  double sum = 0.0;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
    return NAN;

  for (i = 0; i < number_of_variables; ++i) {
    double exponent = 2.0 + (4.0 * (double)(long)i) / ((double)(long)number_of_variables - 1.0);
    sum += pow(fabs(x[i]), exponent);
  }
  result = sqrt(sum);

  return result;
}

/**
 * @brief Implements the sign function.
 */
double sign(double x) {

  if (x > 0)
    return 1;
  if (x < 0)
    return -1;
  return 0;
}

/**
 * @brief The different powers function.
 */
struct f_different_powers_t final : public coco_problem_s {
  explicit f_different_powers_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("different powers function", -5.0, 5.0, 0.0);
    set_id("%s_d%02lu", "different_powers", number_of_variables);

    /* Compute best solution */
    evaluate_function(best_parameter.data(), best_value.data());
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_different_powers_raw(x, number_of_variables);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the gradient of the function "different powers".
   */
  void evaluate_gradient(const double *x, double *y) override {

    size_t i;
    double sum = 0.0;
    double aux;

    for (i = 0; i < number_of_variables; ++i) {
      aux = 2.0 + (4.0 * (double)(long)i) / ((double)(long)number_of_variables - 1.0);
      sum += pow(fabs(x[i]), aux);
    }

    for (i = 0; i < number_of_variables; ++i) {
      aux = 2.0 + (4.0 * (double)(long)i) / ((double)(long)number_of_variables - 1.0);
      y[i] = 0.5 * (aux) / (sum);
      aux -= 1.0;
      y[i] *= pow(fabs(x[i]), aux) * sign(x[i]);
    }
  }
};

/**
 * @brief Allocates the basic different powers problem.
 */
static coco_problem_ptr_t f_different_powers_allocate(const size_t number_of_variables) {
  return std::make_unique<f_different_powers_t>(number_of_variables);
}

/**
 * @brief Creates the BBOB different powers problem.
 */
static coco_problem_ptr_t f_different_powers_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                   const size_t instance, const long rseed,
                                                                   const char *problem_id_template,
                                                                   const char *problem_name_template) {

  double fopt;
  coco_problem_ptr_t problem;

  coco_matrix_t M(dimension, dimension);
  std::vector<double> b(dimension);

  std::vector<double> xopt(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt.data(), rseed, dimension);
  } else {
    bbob2009_compute_xopt(xopt.data(), rseed, dimension);
  }

  const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
  M = rot1;

  problem = f_different_powers_allocate(dimension);
  problem = transform_obj_shift(std::move(problem), fopt);
  problem = transform_vars_affine(std::move(problem), M, b);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("3-ill-conditioned");

  return problem;
}

/**
 * @brief Creates the bbob-constrained different powers problem.
 */
static coco_problem_ptr_t f_different_powers_bbob_constrained_problem_allocate(const size_t function,
                                                                               const size_t dimension,
                                                                               const size_t instance, const long rseed,
                                                                               const char *problem_id_template,
                                                                               const char *problem_name_template) {
  /* Different powers function used in bbob-constrained test suite.
   * In this version, the (unconstrained) optimum, xopt, is set to
   * a distance of 1e-2 to the origin. By doing so, the optimum of
   * the constrained problem is at a "reasonable" distance from
   * the unconstrained one and, hence, the constrained problem is not too easy.
   */

  size_t i;
  double fopt, result;
  coco_problem_ptr_t problem;

  coco_matrix_t M(dimension, dimension);
  std::vector<double> b(dimension);

  std::vector<double> xopt(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt.data(), rseed, dimension);
  } else {
    bbob2009_compute_xopt(xopt.data(), rseed, dimension);
  }

  /* Compute Euclidean norm of xopt */
  /* CAVEAT: this implementation is not ideal for large dimensions */
  result = 0.0;
  for (i = 0; i < dimension; ++i) {
    result += xopt[i] * xopt[i];
  }
  result = sqrt(result);

  /* Scale xopt such that the distance to the origin is 1e-2 */
  for (i = 0; i < dimension; ++i) {
    xopt[i] *= 1e-2 / result;
  }

  const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
  M = rot1;

  problem = f_different_powers_allocate(dimension);
  problem = transform_obj_shift(std::move(problem), fopt);
  problem = transform_vars_affine(std::move(problem), M, b);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("3-ill-conditioned");

  return problem;
}

/**
 * @brief Creates the BBOB generalized permuted block-rotated sum of different
 * powers problem.
 */
static coco_problem_ptr_t f_different_powers_permblockdiag_bbob_problem_allocate(
    const size_t function, const size_t dimension, const size_t instance, const long rseed,
    const char *problem_id_template, const char *problem_name_template) {
  double fopt;
  coco_problem_ptr_t problem;

  std::vector<size_t> P1(dimension);
  std::vector<size_t> P2(dimension);

  size_t nb_blocks;
  size_t swap_range;
  size_t nb_swaps;

  std::vector<size_t> block_sizes = coco_get_block_sizes(&nb_blocks, dimension, "bbob-largescale");
  swap_range = coco_get_swap_range(dimension, "bbob-largescale");
  nb_swaps = coco_get_nb_swaps(dimension, "bbob-largescale");

  std::vector<double> xopt(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt.data(), rseed, dimension);
  } else {
    bbob2009_compute_xopt(xopt.data(), rseed, dimension);
  }

  const coco_block_diagonal_matrix_t B = coco_compute_blockrotation(rseed + 1000000, block_sizes.data(), nb_blocks);

  coco_compute_truncated_uniform_swap_permutation(
      P1.data(), rseed + 2000000, dimension, nb_swaps, swap_range);
  coco_compute_truncated_uniform_swap_permutation(
      P2.data(), rseed + 3000000, dimension, nb_swaps, swap_range);

  problem = f_different_powers_allocate(dimension);
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
