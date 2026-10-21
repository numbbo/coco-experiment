/**
 * @file f_step_ellipsoid.cpp
 * @brief Implementation of the step ellipsoid function and problem.
 *
 * The BBOB step ellipsoid function intertwines the variable and objective
 * transformations in such a way that it is hard to devise a composition of
 * generic transformations to implement it. In the end one would have to
 * implement several custom transformations which would be used solely by this
 * problem. Therefore we opt to implement it as a monolithic function instead.
 *
 * TODO: It would be nice to have a generic step ellipsoid function to
 * complement this one.
 */
#include <assert.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "coco_utilities.cpp"
#include "suite_bbob_legacy_code.cpp"

#include "transform_obj_norm_by_dim.cpp"
#include "transform_obj_penalize.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_blockrotation.cpp"
#include "transform_vars_permutation.cpp"
#include "transform_vars_round_step.cpp"
#include "transform_obj_norm_by_dim.cpp"

/**
 * @brief The BBOB step ellipsoid function.
 */
struct f_step_ellipsoid_t final : public coco_problem_s {
  std::vector<double> z_hat, z;
  std::vector<double> xopt;
  double fopt{}, penalty_scale{};
  coco_matrix_t rot1, rot2;

  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double raw(const double *x);

  explicit f_step_ellipsoid_t(const size_t dimension) : coco_problem_s(dimension, 1, 0) {
    init_from_scalars("step ellipsoid function", -5.0, 5.0, 0);
    /* Allocate temporary storage and space for the rotation matrices */
    z_hat.resize(dimension);
    z.resize(dimension);
    xopt.resize(dimension);
    rot1 = coco_matrix_t(dimension, dimension);
    rot2 = coco_matrix_t(dimension, dimension);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = raw(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Implements the step ellipsoid function without connections to any COCO
 * structures.
 */
double f_step_ellipsoid_t::raw(const double *x) {

  static const double condition = 100;
  static const double alpha = 10.0;
  size_t i, j;
  double penalty = 0.0, x1;
  double result;

  assert(number_of_variables > 1);

  if (coco_vector_contains_nan(x, number_of_variables))
    return NAN;

  for (i = 0; i < number_of_variables; ++i) {
    double tmp;
    tmp = fabs(x[i]) - 5.0;
    if (tmp > 0.0)
      penalty += tmp * tmp;
  }
  penalty = penalty * (penalty_scale);
  for (i = 0; i < number_of_variables; ++i) {
    double c1;
    z_hat[i] = 0.0;
    c1 = sqrt(pow(condition / 10., (double)i / (double)(number_of_variables - 1)));
    for (j = 0; j < number_of_variables; ++j) {
      z_hat[i] += c1 * rot2[i][j] * (x[j] - xopt[j]);
    }
  }
  x1 = z_hat[0];

  for (i = 0; i < number_of_variables; ++i) {
    if (fabs(z_hat[i]) > 0.5) /* TODO: Documentation: no fabs() in documentation */
      z_hat[i] = coco_double_round(z_hat[i]);
    else
      z_hat[i] = coco_double_round(alpha * z_hat[i]) / alpha;
  }

  for (i = 0; i < number_of_variables; ++i) {
    z[i] = 0.0;
    for (j = 0; j < number_of_variables; ++j) {
      z[i] += rot1[i][j] * z_hat[j];
    }
  }

  /* Computation core */
  result = 0.0;
  for (i = 0; i < number_of_variables; ++i) {
    double exponent;
    exponent = (double)(long)i / ((double)(long)number_of_variables - 1.0);
    result += pow(condition, exponent) * z[i] * z[i];
  }
  result = 0.1 * coco_double_max(fabs(x1) * 1.0e-4, result) + penalty + fopt;

  return result;
}

/**
 * @brief Creates the BBOB step ellipsoid problem.
 *
 * @note There is no separate basic allocate function.
 */
static coco_problem_ptr_t f_step_ellipsoid_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                 const size_t instance, const long rseed,
                                                                 const double penalty_scale,
                                                                 const char *problem_id_template,
                                                                 const char *problem_name_template) {

  size_t i;
  std::unique_ptr<f_step_ellipsoid_t> problem = std::make_unique<f_step_ellipsoid_t>(dimension);

  problem->fopt = bbob2009_compute_fopt(function, instance);
  problem->penalty_scale = penalty_scale;
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(problem->xopt.data(), rseed, dimension);
  } else {
    bbob2009_compute_xopt(problem->xopt.data(), rseed, dimension);
  }
  bbob2009_compute_rotation(problem->rot1.rows(), rseed + 1000000, dimension);
  bbob2009_compute_rotation(problem->rot2.rows(), rseed, dimension);

  /* Compute best solution
   *
   * OME: Dirty hack for now because I did not want to invert the
   * transformations to find the best_parameter :/
   */
  for (i = 0; i < problem->number_of_variables; i++) {
    problem->best_parameter[i] = problem->xopt[i];
  }
  problem->best_value[0] = problem->fopt;
  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("2-moderate");

  return problem;
}

/**
 * @brief Implements the step ellipsoid function without connections to any COCO
 * structures.
 */
static double f_step_ellipsoid_core(const double *x, const size_t number_of_variables,
                                    f_step_ellipsoid_shared_data_t *shared) {

  static const double condition = 100;
  size_t i;
  double result;
  result = 0.0;

  for (i = 0; i < number_of_variables; ++i) {
    double exponent;
    exponent = (double)(long)i / ((double)(long)number_of_variables - 1.0);
    result += pow(condition, exponent) * x[i] * x[i];
  }
  result = 0.1 * coco_double_max(
                     shared->zhat_1 * 1.0e-4, result);
  return result;
}

/**
 * @brief The step ellipsoid function of the large-scale suite.
 *
 * Shares \hat{z}_1 with transform_vars_round_step, which sets it.
 */
struct f_step_ellipsoid_permblock_t final : public coco_problem_s {
  f_step_ellipsoid_shared_data_t shared;

  explicit f_step_ellipsoid_permblock_t(const size_t number_of_variables)
      : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("step ellipsoid function", -5.0, 5.0, 0.0);
    shared.zhat_1 = 0; /*needed for xopt evaluation*/
    set_id("%s_d%02lu", "step_ellipsoid", number_of_variables);

    /* Compute best solution, here done outside after the zhat is set to the best_value */
    evaluate_function(best_parameter.data(), best_value.data());
  }

  /**
   * @brief Uses the raw function to evaluate the ls COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_step_ellipsoid_core(x, number_of_variables, &shared);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Allocates the basic step ellipsoid problem.
 * an additional coordinate is added that will contain the value of \hat{z}_1 but that is ignored by functions other
 * that f_step_ellipsoid_core and transform_vars_round_step. The latter sets it.
 */
static std::unique_ptr<f_step_ellipsoid_permblock_t> f_step_ellipsoid_allocate(const size_t number_of_variables) {
  return std::make_unique<f_step_ellipsoid_permblock_t>(number_of_variables);
}

/**
 * @brief Creates the BBOB permuted block-rotated step ellipsoid problem.
 *
 * Wassim: TODO: consider implementing it sub-problem style
 * Wassim: TODO: make the zhat1 value default to x1 when no transformation is
 * applied and the data type defined here
 */
static coco_problem_ptr_t f_step_ellipsoid_permblockdiag_bbob_problem_allocate(const size_t function,
                                                                               const size_t dimension,
                                                                               const size_t instance, const long rseed,
                                                                               const char *problem_id_template,
                                                                               const char *problem_name_template) {
  double alpha = 10.; /*parameter of rounding*/
  double fopt;
  coco_problem_ptr_t problem;
  size_t nb_blocks1, nb_blocks2, swap_range1, swap_range2, nb_swaps1, nb_swaps2;
  double penalty_factor = 1.;

  std::vector<double> xopt(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt.data(), rseed, dimension);
  } else {
    bbob2009_compute_xopt(xopt.data(), rseed, dimension);
  }

  std::vector<size_t> block_sizes1 = coco_get_block_sizes(&nb_blocks1, dimension, "bbob-largescale");
  std::vector<size_t> block_sizes2 = coco_get_block_sizes(&nb_blocks2, dimension, "bbob-largescale");
  swap_range1 = coco_get_swap_range(dimension, "bbob-largescale");
  swap_range2 = coco_get_swap_range(dimension, "bbob-largescale");
  nb_swaps1 = coco_get_nb_swaps(dimension, "bbob-largescale");
  nb_swaps2 = coco_get_nb_swaps(dimension, "bbob-largescale");

  coco_matrix_t B1 = coco_block_matrix(dimension, block_sizes1.data(), nb_blocks1);
  coco_matrix_t B2 = coco_block_matrix(dimension, block_sizes2.data(), nb_blocks2);
  coco_compute_blockrotation(B1.rows(), rseed + 1000000, dimension, block_sizes1.data(),
                             nb_blocks1);
  coco_compute_blockrotation(B2.rows(), rseed, dimension, block_sizes2.data(), nb_blocks2);

  std::vector<size_t> P11(dimension);
  std::vector<size_t> P12(dimension);
  std::vector<size_t> P21(dimension);
  std::vector<size_t> P22(dimension);
  coco_compute_truncated_uniform_swap_permutation(
      P11.data(), rseed + 2000000, dimension, nb_swaps1, swap_range1);
  coco_compute_truncated_uniform_swap_permutation(
      P12.data(), rseed + 3000000, dimension, nb_swaps1, swap_range1);
  coco_compute_truncated_uniform_swap_permutation(
      P21.data(), rseed + 4000000, dimension, nb_swaps2, swap_range2);
  coco_compute_truncated_uniform_swap_permutation(
      P22.data(), rseed + 5000000, dimension, nb_swaps2, swap_range2);

  std::unique_ptr<f_step_ellipsoid_permblock_t> step_ellipsoid = f_step_ellipsoid_allocate(dimension);
  f_step_ellipsoid_permblock_t *raw_problem = step_ellipsoid.get();
  problem = std::move(step_ellipsoid);

  problem = transform_vars_permutation(std::move(problem), P22.data(), dimension);
  problem = transform_vars_blockrotation(std::move(problem), B1.rows(), dimension, block_sizes1.data(), nb_blocks1);
  problem = transform_vars_permutation(std::move(problem), P21.data(), dimension);
  problem = transform_vars_round_step(std::move(problem), alpha, &raw_problem->shared);

  problem = transform_vars_conditioning(std::move(problem), 10.0);
  problem = transform_vars_permutation(std::move(problem), P12.data(), dimension);
  problem = transform_vars_blockrotation(std::move(problem), B2.rows(), dimension, block_sizes2.data(), nb_blocks2);
  problem = transform_vars_permutation(std::move(problem), P11.data(), dimension);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);

  problem = transform_obj_norm_by_dim(std::move(problem));
  problem = transform_obj_penalize(std::move(problem), penalty_factor);
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("2-moderate");

  return problem;
}
