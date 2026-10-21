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
static double f_katsuura_raw(const double *x,
                             const size_t number_of_variables) {

  size_t i, j;
  double tmp, tmp2;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
    return NAN;

  /* Computation core */
  result = 1.0;
  for (i = 0; i < number_of_variables; ++i) {
    tmp = 0;
    for (j = 1; j < 33; ++j) {
      tmp2 = pow(2., (double)j);
      tmp += fabs(tmp2 * x[i] - coco_double_round(tmp2 * x[i])) / tmp2;
    }
    tmp = 1.0 + ((double)(long)i + 1) * tmp;
    /*result *= tmp;*/ /* Wassim TODO: delete once consistency check passed*/
    result *= pow(tmp, 10. / pow((double)number_of_variables, 1.2));
  }
  /*result = 10. / ((double) number_of_variables) / ((double) number_of_variables)
   * (-1. + pow(result, 10. / pow((double) number_of_variables, 1.2)));*/
  result = 10. / ((double)number_of_variables) / ((double)number_of_variables) * (-1. + result);

  return result;
}

/**
 * @brief The Katsuura function.
 */
struct f_katsuura_t final : public coco_problem_s {
  explicit f_katsuura_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    coco_problem_init_from_scalars(this, "Katsuura function", -5.0, 5.0, 0.0);
    coco_problem_set_id(this, "%s_d%02lu", "katsuura", number_of_variables);

    /* Compute best solution */
    evaluate_function(best_parameter.data(), best_value.data());
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_katsuura_raw(x, number_of_variables);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Allocates the basic Katsuura problem.
 */
static coco_problem_t *f_katsuura_allocate(const size_t number_of_variables) {
  return new f_katsuura_t(number_of_variables);
}

/**
 * @brief Creates the BBOB Katsuura problem.
 */
static coco_problem_t *f_katsuura_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                        const size_t instance, const long rseed,
                                                        const char *problem_id_template,
                                                        const char *problem_name_template) {
  double fopt;
  coco_problem_t *problem = nullptr;
  size_t i, j, k;
  std::vector<double> M(dimension * dimension);
  std::vector<double> b(dimension);
  double *current_row;

  const double penalty_factor = 1.0;

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
        current_row[j] += rot1[i][k] * pow(sqrt(100), exponent) * rot2[k][j];
      }
    }
  }

  problem = f_katsuura_allocate(dimension);
  problem = transform_obj_shift(
      problem, fopt); /*There is no shift 'fopt' in the definition*/
  problem = transform_vars_affine(problem, M.data(), b.data(), dimension);
  problem = transform_vars_shift(problem, xopt.data(), 0);
  problem = transform_obj_penalize(problem, penalty_factor);


  coco_problem_set_id(problem, problem_id_template, function, instance,
                      dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance,
                        dimension);
  coco_problem_set_type(problem, "5-weakly-structured");

  return problem;
}

/**
 * @brief Creates the BBOB permuted block-rotated Katsuura problem.
 */
static coco_problem_t *f_katsuura_permblockdiag_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                      const size_t instance, const long rseed,
                                                                      const char *problem_id_template,
                                                                      const char *problem_name_template) {
  double fopt;
  coco_problem_t *problem = nullptr;
  const double penalty_factor = 1.0;
  std::vector<size_t> P11(dimension);
  std::vector<size_t> P12(dimension);
  std::vector<size_t> P21(dimension);
  std::vector<size_t> P22(dimension);
  size_t nb_blocks1;
  size_t nb_blocks2;
  size_t swap_range1;
  size_t swap_range2;
  size_t nb_swaps1;
  size_t nb_swaps2;

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

  coco_compute_truncated_uniform_swap_permutation(
      P11.data(), rseed + 3000000, dimension, nb_swaps1, swap_range1);
  coco_compute_truncated_uniform_swap_permutation(
      P12.data(), rseed + 4000000, dimension, nb_swaps1, swap_range1);
  coco_compute_truncated_uniform_swap_permutation(
      P21.data(), rseed + 5000000, dimension, nb_swaps2, swap_range2);
  coco_compute_truncated_uniform_swap_permutation(
      P22.data(), rseed + 6000000, dimension, nb_swaps2, swap_range2);

  problem = f_katsuura_allocate(dimension);
  problem = transform_vars_permutation(problem, P22.data(), dimension);
  problem = transform_vars_blockrotation(problem, B1.rows(), dimension,
                                         block_sizes2.data(), nb_blocks2);
  problem = transform_vars_permutation(problem, P21.data(), dimension);
  problem = transform_vars_conditioning(problem, 100.0);
  problem = transform_vars_permutation(problem, P12.data(), dimension);
  problem = transform_vars_blockrotation(problem, B2.rows(), dimension,
                                         block_sizes1.data(), nb_blocks1);
  problem = transform_vars_permutation(problem, P11.data(), dimension);
  problem = transform_vars_shift(problem, xopt.data(), 0);

  /*problem = transform_obj_norm_by_dim(problem);*/ /* Wassim: does not seem to
                                                       be needed*/
  problem = transform_obj_penalize(problem, penalty_factor);
  problem =
      transform_obj_shift(problem, fopt); /*TODO: documentation, there is no fopt in the definition of this function*/

  coco_problem_set_id(problem, problem_id_template, function, instance,
                      dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance,
                        dimension);
  coco_problem_set_type(problem, "5-weakly-structured");

  return problem;
}
