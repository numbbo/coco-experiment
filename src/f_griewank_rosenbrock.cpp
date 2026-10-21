/**
 * @file f_griewank_rosenbrock.cpp
 * @brief Implementation of the Griewank-Rosenbrock function and problem.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_norm_by_dim.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_affine.cpp"
#include "transform_vars_blockrotation.cpp"
#include "transform_vars_permutation.cpp"
#include "transform_vars_scale.cpp"
#include "transform_vars_shift.cpp"

/**
 * @brief Data type for the griewank rosenbrock problem
 */
struct f_griewank_rosenbrock_data_t {
  double facftrue;
};

/**
 * @brief Implements the Griewank-Rosenbrock function without connections to any
 * COCO structures.
 */
static double f_griewank_rosenbrock_raw(const double *x, const size_t number_of_variables,
                                        f_griewank_rosenbrock_data_t *data) {

  size_t i = 0;
  double tmp = 0;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
    return NAN;

  /* Computation core */
  result = 0.0;
  for (i = 0; i < number_of_variables - 1; ++i) {
    const double c1 = x[i] * x[i] - x[i + 1];
    const double c2 = 1.0 - x[i];
    tmp = 100.0 * c1 * c1 + c2 * c2;
    result += tmp / 4000. - cos(tmp);
  }
  result = data->facftrue + data->facftrue * result / (double)(number_of_variables - 1);

  return result;
}

/**
 * @brief The Griewank-Rosenbrock function.
 */
struct f_griewank_rosenbrock_t final : public coco_problem_s {
  f_griewank_rosenbrock_data_t params;

  f_griewank_rosenbrock_t(const size_t number_of_variables, double facftrue)
      : coco_problem_s(number_of_variables, 1, 0) {
    coco_problem_init_from_scalars(this, "Griewank Rosenbrock function", -5.0, 5.0, 1);
    coco_problem_set_id(this, "%s_d%02lu", "griewank_rosenbrock", number_of_variables);

    params.facftrue = facftrue;

    /* Compute best solution */
    evaluate_function(best_parameter.data(), best_value.data());
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_griewank_rosenbrock_raw(x, number_of_variables, &params);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Allocates the basic Griewank-Rosenbrock problem.
 */
static coco_problem_t *f_griewank_rosenbrock_allocate(const size_t number_of_variables, double facftrue) {
  return new f_griewank_rosenbrock_t(number_of_variables, facftrue);
}

/**
 * @brief Creates the BBOB Griewank-Rosenbrock problem.
 */
static coco_problem_t *f_griewank_rosenbrock_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                   const size_t instance, const long rseed,
                                                                   const void *args, const char *problem_id_template,
                                                                   const char *problem_name_template) {
  double fopt;
  coco_problem_t *problem = nullptr;
  size_t i, j;

  std::vector<double> M(dimension * dimension);
  std::vector<double> b(dimension);
  std::vector<double> shift(dimension);
  double scales;
  double tmp; /* Wassim: will serve to set the optimal solution "manually"*/

  fopt = bbob2009_compute_fopt(function, instance);
  for (i = 0; i < dimension; ++i) {
    shift[i] = -0.5;
  }

  coco_matrix_t rot1(dimension, dimension);
  bbob2009_compute_rotation(rot1.rows(), rseed, dimension);
  scales = coco_double_max(1., sqrt((double)dimension) / 8.);
  for (i = 0; i < dimension; ++i) {
    for (j = 0; j < dimension; ++j) {
      rot1[i][j] *= scales;
    }
  }

  f_griewank_rosenbrock_args_t *f_griewank_rosenbrock_args;
  f_griewank_rosenbrock_args = ((f_griewank_rosenbrock_args_t *)args);
  problem = f_griewank_rosenbrock_allocate(dimension, f_griewank_rosenbrock_args->facftrue);
  problem = transform_obj_shift(problem, fopt);
  problem = transform_vars_shift(problem, shift.data(), 0);
  bbob2009_copy_rotation_matrix(rot1.rows(), M.data(), b.data(), dimension);

  for (i = 0; i < dimension; i++) {
    problem->best_parameter[i] =
        0; /* Wassim: TODO: not a proper way of avoiding to trigger coco_warning("transform_vars_affine():
              'best_parameter' not updated, set to NAN")*/
  }
  problem = transform_vars_affine(problem, M.data(), b.data(), dimension);
  for (j = 0; j < dimension; ++j) { /* Wassim: manually set xopt = rot1^T
                                       ones(dimension)/(2*factor) */
    tmp = 0;
    for (i = 0; i < dimension; ++i) {
      tmp += rot1[i][j];
    }
    problem->best_parameter[j] = tmp / (2. * scales);
  }

  coco_problem_set_id(problem, problem_id_template, function, instance,
                      dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance,
                        dimension);
  coco_problem_set_type(problem, "4-multi-modal");

  return problem;
}

/**
 * @brief Creates the BBOB permuted block-rotated Griewank-Rosenbrock problem.
 */
static coco_problem_t *f_griewank_rosenbrock_permblockdiag_bbob_bbob_problem_allocate(
    const size_t function, const size_t dimension, const size_t instance, const long rseed, const void *args,
    const char *problem_id_template, const char *problem_name_template) {
  double fopt;
  coco_problem_t *problem = nullptr;
  double scales;
  size_t i, j, k, next_bs_change;

  std::vector<size_t> P1(dimension);
  std::vector<size_t> P2(dimension);
  size_t block_size;
  size_t nb_blocks;
  size_t swap_range;
  size_t nb_swaps;
  double tmp;                                               /* Manh: will serve to set the optimal solution "manually"*/
  std::vector<double> best_parameter(dimension); /* Manh: will serve to set the optimal solution "manually"*/

  std::vector<size_t> block_sizes = coco_get_block_sizes(&nb_blocks, dimension, "bbob-largescale");
  block_size = coco_rotation_matrix_block_size(dimension);
  swap_range = coco_get_swap_range(dimension, "bbob-largescale");
  nb_swaps = coco_get_nb_swaps(dimension, "bbob-largescale");

  fopt = bbob2009_compute_fopt(function, instance);
  scales = coco_double_max(1.0, sqrt((double)block_size) / 8.0);
  std::vector<double> shift(dimension);
  for (i = 0; i < dimension; ++i) {
    shift[i] = -0.5;
  }

  coco_matrix_t B = coco_block_matrix(dimension, block_sizes.data(), nb_blocks);

  coco_compute_blockrotation(B.rows(), rseed, dimension, block_sizes.data(), nb_blocks);
  coco_compute_truncated_uniform_swap_permutation(P1.data(), rseed + 2000000, dimension, nb_swaps, swap_range);
  coco_compute_truncated_uniform_swap_permutation(P2.data(), rseed + 3000000, dimension, nb_swaps, swap_range);

  f_griewank_rosenbrock_args_t *f_griewank_rosenbrock_args;
  f_griewank_rosenbrock_args = ((f_griewank_rosenbrock_args_t *)args);

  problem = f_griewank_rosenbrock_allocate(dimension, f_griewank_rosenbrock_args->facftrue);
  problem = transform_vars_shift(problem, shift.data(), 0);
  problem = transform_vars_scale(problem, scales);
  problem = transform_vars_permutation(problem, P2.data(), dimension);
  problem = transform_vars_blockrotation(problem, B.rows(), dimension,
                                         block_sizes.data(), nb_blocks);
  problem = transform_vars_permutation(problem, P1.data(), dimension);
  
  /*problem = transform_obj_norm_by_dim(problem);*/ /* Wassim: there is already a normalization by dimension*/
  problem = transform_obj_shift(problem, fopt);

  /* Manh: manually set xopt = rot1^T ones(dimension)/(2*scales) */
  next_bs_change = 0;
  for (k = 0; k < nb_blocks; ++k) {
    for (j = 0; j < block_sizes[k]; ++j) { /* Manh: firstly, set xopt_1 = (B^T)*(P_2^T)*ones(dimension)/(2*scales) */
      tmp = 0;
      for (i = 0; i < block_sizes[k]; ++i) {
        tmp += B[next_bs_change + i][j];
      }
      best_parameter[next_bs_change + j] = tmp / (2. * scales);
    }
    next_bs_change += block_sizes[k];
  }

  for (j = 0; j < dimension;
       ++j) { /* Manh: secondly, set xopt = (P_1^T)* xopt_1 */
    problem->best_parameter[P1[j]] = best_parameter[j];
  }

  coco_problem_set_id(problem, problem_id_template, function, instance,
                      dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance,
                        dimension);
  coco_problem_set_type(problem, "4-multi-modal");

  return problem;
}
