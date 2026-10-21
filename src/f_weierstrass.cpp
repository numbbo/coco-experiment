/**
 * @file f_weierstrass.cpp
 * @brief Implementation of the Weierstrass function and problem.
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
#include "transform_vars_blockrotation.cpp"
#include "transform_vars_oscillate.cpp"
#include "transform_vars_permutation.cpp"
#include "transform_vars_shift.cpp"

/** @brief Number of summands in the Weierstrass problem. */
#define F_WEIERSTRASS_SUMMANDS 12

/**
 * @brief The Weierstrass function.
 */
struct f_weierstrass_t final : public coco_problem_s {
  double f0{};
  double ak[F_WEIERSTRASS_SUMMANDS]{};
  double bk[F_WEIERSTRASS_SUMMANDS]{};

  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double raw(const double *x) const;

  explicit f_weierstrass_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    size_t i;

    init_from_scalars("Weierstrass function", -5.0, 5.0, 0.);
    set_id("%s_d%02lu", "weierstrass", number_of_variables);

    f0 = 0.0;
    for (i = 0; i < F_WEIERSTRASS_SUMMANDS; ++i) {
      ak[i] = pow(0.5, (double)i);
      bk[i] = pow(3., (double)i);
      f0 += ak[i] * cos(2 * coco_pi * bk[i] * 0.5);
    }

    /* Compute best solution */
    std::vector<double> non_unique_best_value(number_of_variables, 0.0);
    evaluate_function(non_unique_best_value.data(), best_value.data());
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
 * @brief Implements the Weierstrass function without connections to any COCO
 * structures.
 */
double f_weierstrass_t::raw(const double *x) const {

  size_t i, j;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
    return NAN;

  result = 0.0;
  for (i = 0; i < number_of_variables; ++i) {
    for (j = 0; j < F_WEIERSTRASS_SUMMANDS; ++j) {
      result += cos(2 * coco_pi * (x[i] + 0.5) * bk[j]) * ak[j];
    }
  }
  result = 10.0 * pow(result / (double)(long)number_of_variables - f0, 3.0);

  return result;
}

/**
 * @brief Allocates the basic Weierstrass problem.
 */
static coco_problem_ptr_t f_weierstrass_allocate(const size_t number_of_variables) {
  return std::make_unique<f_weierstrass_t>(number_of_variables);
}

/**
 * @brief Creates the BBOB Weierstrass problem.
 */
static coco_problem_ptr_t f_weierstrass_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                              const size_t instance, const long rseed,
                                                              const char *problem_id_template,
                                                              const char *problem_name_template) {
  double fopt;
  coco_problem_ptr_t problem;
  size_t i, j, k;
  std::vector<double> M(dimension * dimension);
  std::vector<double> b(dimension);
  double *current_row;

  const double condition = 100.0;
  const double penalty_factor = 10.0 / (double)dimension;

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
        const double base = 1.0 / sqrt(condition);
        const double exponent = 1.0 * (int)k / ((double)(long)dimension - 1.0);
        current_row[j] += rot1[i][k] * pow(base, exponent) * rot2[k][j];
      }
    }
  }

  problem = f_weierstrass_allocate(dimension);
  problem = transform_obj_shift(std::move(problem), fopt);
  problem = transform_vars_affine(std::move(problem), M.data(), b.data(), dimension);
  problem = transform_vars_oscillate(std::move(problem));
  bbob2009_copy_rotation_matrix(rot1.rows(), M.data(), b.data(), dimension);
  problem = transform_vars_affine(std::move(problem), M.data(), b.data(), dimension);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);
  problem = transform_obj_penalize(std::move(problem), penalty_factor);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("4-multi-modal");

  return problem;
}

/**
 * @brief Creates the BBOB permuted block-rotated Weierstrass problem.
 */
static coco_problem_ptr_t f_weierstrass_permblockdiag_bbob_problem_allocate(const size_t function,
                                                                            const size_t dimension,
                                                                            const size_t instance, const long rseed,
                                                                            const char *problem_id_template,
                                                                            const char *problem_name_template) {
  double fopt;
  coco_problem_ptr_t problem;
  std::vector<size_t> P11(dimension);
  std::vector<size_t> P21(dimension);
  std::vector<size_t> P12(dimension);
  std::vector<size_t> P22(dimension);
  size_t nb_blocks1;
  size_t nb_blocks2;
  size_t swap_range;
  size_t nb_swaps;
  const double condition = 100.0, penalty_factor = 10.0 / (double)dimension;

  std::vector<size_t> block_sizes1 = coco_get_block_sizes(&nb_blocks1, dimension, "bbob-largescale");
  std::vector<size_t> block_sizes2 = coco_get_block_sizes(&nb_blocks2, dimension, "bbob-largescale");
  swap_range = coco_get_swap_range(dimension, "bbob-largescale");
  nb_swaps = coco_get_nb_swaps(dimension, "bbob-largescale");

  std::vector<double> xopt(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt.data(), rseed, dimension);
  } else {
    bbob2009_compute_xopt(xopt.data(), rseed, dimension);
  }

  coco_matrix_t B1 = coco_block_matrix(dimension, block_sizes1.data(), nb_blocks1);
  coco_matrix_t B2 = coco_block_matrix(dimension, block_sizes2.data(), nb_blocks2);

  coco_compute_blockrotation(B1.rows(), rseed + 1000000, dimension, block_sizes1.data(), nb_blocks1);
  coco_compute_blockrotation(B2.rows(), rseed, dimension, block_sizes2.data(), nb_blocks2);

  coco_compute_truncated_uniform_swap_permutation(P11.data(), rseed + 3000000, dimension, nb_swaps, swap_range);
  coco_compute_truncated_uniform_swap_permutation(P21.data(), rseed + 4000000, dimension, nb_swaps, swap_range);
  coco_compute_truncated_uniform_swap_permutation(P12.data(), rseed + 5000000, dimension, nb_swaps, swap_range);
  coco_compute_truncated_uniform_swap_permutation(P22.data(), rseed + 6000000, dimension, nb_swaps, swap_range);

  problem = f_weierstrass_allocate(dimension);
  problem = transform_vars_permutation(std::move(problem), P22.data(), dimension);
  problem = transform_vars_blockrotation(std::move(problem), B1.rows(), dimension, block_sizes1.data(), nb_blocks1);
  problem = transform_vars_permutation(std::move(problem), P12.data(), dimension);

  problem = transform_vars_conditioning(std::move(problem), 1.0 / condition);
  problem = transform_vars_permutation(std::move(problem), P21.data(), dimension);
  problem = transform_vars_blockrotation(std::move(problem), B2.rows(), dimension, block_sizes2.data(), nb_blocks2);
  problem = transform_vars_permutation(std::move(problem), P11.data(), dimension);

  problem = transform_vars_oscillate(std::move(problem));
  problem = transform_vars_permutation(std::move(problem), P22.data(), dimension);
  problem = transform_vars_blockrotation(std::move(problem), B1.rows(), dimension, block_sizes1.data(), nb_blocks1);
  problem = transform_vars_permutation(std::move(problem), P12.data(), dimension);

  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);
  /*problem = transform_obj_norm_by_dim(std::move(problem));*/ /* Wassim: there is already a normalization by
                                                                  dimension*/
  problem = transform_obj_penalize(std::move(problem), penalty_factor);
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("4-multi-modal");

  return problem;
}

#undef F_WEIERSTRASS_SUMMANDS
