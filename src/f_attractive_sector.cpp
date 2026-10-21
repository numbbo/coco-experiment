/**
 * @file f_attractive_sector.cpp
 * @brief Implementation of the attractive sector function and problem.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_norm_by_dim.cpp"
#include "transform_obj_oscillate.cpp"
#include "transform_obj_power.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_affine.cpp"
#include "transform_vars_blockrotation.cpp"
#include "transform_vars_conditioning.cpp"
#include "transform_vars_permutation.cpp"
#include "transform_vars_shift.cpp"

/**
 * @brief The attractive sector function.
 */
struct f_attractive_sector_t final : public coco_problem_s {
  std::vector<double> xopt;

  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double raw(const double *x) const;

  f_attractive_sector_t(const size_t number_of_variables, const double *xopt)
      : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("attractive sector function", -5.0, 5.0, 0.0);
    set_id("%s_d%02lu", "attractive_sector", number_of_variables);

    this->xopt.assign(xopt, xopt + number_of_variables);

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
};

/**
 * @brief Implements the attractive sector function without connections to any
 * COCO structures.
 */
double f_attractive_sector_t::raw(const double *x) const {
  size_t i;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
    return NAN;

  result = 0.0;
  for (i = 0; i < number_of_variables; ++i) {
    if (xopt[i] * x[i] > 0.0) {
      result += 100.0 * 100.0 * x[i] * x[i];
    } else {
      result += x[i] * x[i];
    }
  }
  return result;
}

/**
 * @brief Allocates the basic attractive sector problem.
 */
static coco_problem_ptr_t f_attractive_sector_allocate(const size_t number_of_variables, const double *xopt) {
  return std::make_unique<f_attractive_sector_t>(number_of_variables, xopt);
}

/**
 * @brief Creates the BBOB attractive sector problem.
 */
static coco_problem_ptr_t f_attractive_sector_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                    const size_t instance, const long rseed,
                                                                    const char *problem_id_template,
                                                                    const char *problem_name_template) {
  double fopt;
  coco_problem_ptr_t problem;
  size_t k;
  coco_matrix_t M(dimension, dimension);
  std::vector<double> b(dimension);

  std::vector<double> xopt(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt.data(), rseed, dimension);
  } else {
    bbob2009_compute_xopt(xopt.data(), rseed, dimension);
  }

  /* Compute affine transformation M from two rotation matrices */
  const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
  const coco_matrix_t rot2 = bbob2009_compute_rotation(rseed, dimension);
  std::vector<double> scales(dimension);
  for (k = 0; k < dimension; ++k) {
    const double exponent = 1.0 * (int)k / ((double)(long)dimension - 1.0);
    scales[k] = pow(sqrt(10.0), exponent);
  }
  M = coco_matrix_product_with_diagonal(rot1, scales, rot2);

  problem = f_attractive_sector_allocate(dimension, xopt.data());
  problem = transform_obj_oscillate(std::move(problem));
  problem = transform_obj_power(std::move(problem), 0.9);
  problem = transform_obj_shift(std::move(problem), fopt);
  problem = transform_vars_affine(std::move(problem), M, b);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("2-moderate");

  return problem;
}

static coco_problem_ptr_t f_attractive_sector_permblockdiag_bbob_problem_allocate(
    const size_t function, const size_t dimension, const size_t instance, const long rseed,
    const char *problem_id_template, const char *problem_name_template) {

  double fopt;
  coco_problem_ptr_t problem;
  const double condition = 10.0;
  /* each of R and Q might have its own parameter values */
  size_t nb_blocks1, nb_blocks2;
  size_t swap_range1, swap_range2;
  size_t nb_swaps1, nb_swaps2;

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

  const coco_block_diagonal_matrix_t B1 = coco_compute_blockrotation(rseed + 1000000, block_sizes1.data(), nb_blocks1);
  const coco_block_diagonal_matrix_t B2 = coco_compute_blockrotation(rseed, block_sizes2.data(), nb_blocks2);

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

  problem = f_attractive_sector_allocate(dimension, xopt.data());
  problem = transform_obj_norm_by_dim(std::move(problem));
  problem = transform_obj_oscillate(std::move(problem));
  problem = transform_obj_power(std::move(problem), 0.9);
  problem = transform_obj_shift(std::move(problem), fopt);

  problem = transform_vars_permutation(std::move(problem), P22.data(), dimension);
  problem = transform_vars_blockrotation(std::move(problem), B1);
  problem = transform_vars_permutation(std::move(problem), P21.data(), dimension);
  problem = transform_vars_conditioning(std::move(problem), condition);
  problem = transform_vars_permutation(std::move(problem), P12.data(), dimension);
  problem = transform_vars_blockrotation(std::move(problem), B2);
  problem = transform_vars_permutation(std::move(problem), P11.data(), dimension);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("2-moderate");

  return problem;
}
