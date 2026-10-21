/**
 * @file f_ellipsoid.cpp
 * @brief Implementation of the ellipsoid function and problem.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "transform_obj_norm_by_dim.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_affine.cpp"
#include "transform_vars_blockrotation.cpp"
#include "transform_vars_oscillate.cpp"
#include "transform_vars_permutation.cpp"
#include "transform_vars_blockrotation.cpp"
#include "transform_obj_norm_by_dim.cpp"

/**
 * @brief Data type for the ellipsoid problem.
 */
struct f_ellipsoid_data_t {
  double conditioning;
};

/**
 * @brief Implements the ellipsoid function without connections to any COCO
 * structures.
 */
static double f_ellipsoid_raw(const double *x, const size_t number_of_variables, f_ellipsoid_data_t *data) {

  size_t i = 0;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
    return NAN;

  result = x[i] * x[i];
  for (i = 1; i < number_of_variables; ++i) {
    const double exponent = 1.0 * (double)(long)i / ((double)(long)number_of_variables - 1.0);
    result += pow(data->conditioning, exponent) * x[i] * x[i];
  }

  return result;
}

/**
 * @brief The ellipsoid function.
 */
struct f_ellipsoid_t final : public coco_problem_s {
  f_ellipsoid_data_t params;

  f_ellipsoid_t(const size_t number_of_variables, const double conditioning)
      : coco_problem_s(number_of_variables, 1, 0) {
    coco_problem_init_from_scalars(this, "ellipsoid function", -5.0, 5.0, 0.0);
    coco_problem_set_id(this, "%s_d%02lu", "ellipsoid", number_of_variables);

    params.conditioning = conditioning;

    /* Compute best solution */
    evaluate_function(best_parameter.data(), best_value.data());
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_ellipsoid_raw(x, number_of_variables, &params);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the gradient of the ellipsoid function.
   */
  void evaluate_gradient(const double *x, double *y) override {
    double exponent;
    size_t i = 0;

    for (i = 0; i < number_of_variables; ++i) {
      exponent = 1.0 * (double)(long)i / ((double)(long)number_of_variables - 1.0);
      y[i] = 2.0 * pow(params.conditioning, exponent) * x[i];
    }
  }
};

/**
 * @brief Allocates the basic ellipsoid problem.
 */
static coco_problem_ptr_t f_ellipsoid_allocate(const size_t number_of_variables, const double conditioning) {
  return std::make_unique<f_ellipsoid_t>(number_of_variables, conditioning);
}

/**
 * @brief Creates the BBOB ellipsoid problem.
 */
static coco_problem_ptr_t f_ellipsoid_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                            const size_t instance, const long rseed,
                                                            const char *problem_id_template,
                                                            const char *problem_name_template) {
  double fopt;
  coco_problem_ptr_t problem;

  std::vector<double> xopt(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt.data(), rseed, dimension);
  } else {
    bbob2009_compute_xopt(xopt.data(), rseed, dimension);
  }

  problem = f_ellipsoid_allocate(dimension, 1.0e6);
  problem = transform_vars_oscillate(std::move(problem));
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);

  /*if large scale test-bed, normalize by dim*/
  if (coco_strfind(problem_name_template, "BBOB large-scale suite") >= 0) {
    problem = transform_obj_norm_by_dim(std::move(problem));
  }
  problem = transform_obj_shift(std::move(problem), fopt);

  coco_problem_set_id(problem.get(), problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem.get(), problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem.get(), "1-separable");

  return problem;
}

/**
 * @brief Creates the BBOB rotated ellipsoid problem.
 */
static coco_problem_ptr_t f_ellipsoid_rotated_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                    const size_t instance, const long rseed,
                                                                    const void *args, const char *problem_id_template,
                                                                    const char *problem_name_template) {
  double fopt;
  coco_problem_ptr_t problem;

  f_ellipsoid_args_t *f_ellipsoid_args;
  f_ellipsoid_args = ((f_ellipsoid_args_t *)args);
  std::vector<double> M(dimension * dimension);
  std::vector<double> b(dimension);

  double conditioning = f_ellipsoid_args->conditioning;
  std::vector<double> xopt(dimension);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt.data(), rseed, dimension);
  } else {
    bbob2009_compute_xopt(xopt.data(), rseed, dimension);
  }
  fopt = bbob2009_compute_fopt(function, instance);

  coco_matrix_t rot1(dimension, dimension);
  bbob2009_compute_rotation(rot1.rows(), rseed + 1000000, dimension);
  bbob2009_copy_rotation_matrix(rot1.rows(), M.data(), b.data(), dimension);
  problem = f_ellipsoid_allocate(dimension, conditioning);
  problem = transform_vars_oscillate(std::move(problem));
  problem = transform_vars_affine(std::move(problem), M.data(), b.data(), dimension);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);
  problem = transform_obj_shift(std::move(problem), fopt);
  coco_problem_set_id(problem.get(), problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem.get(), problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem.get(), "3-ill-conditioned");

  return problem;
}

/**
 * @brief Creates the BBOB permuted block-rotated ellipsoid problem.
 */
static coco_problem_ptr_t f_ellipsoid_permblockdiag_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                          const size_t instance, const long rseed,
                                                                          const char *problem_id_template,
                                                                          const char *problem_name_template) {
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
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt.data(), rseed, dimension);
  } else {
    bbob2009_compute_xopt(xopt.data(), rseed, dimension);
  }
  fopt = bbob2009_compute_fopt(function, instance);

  coco_matrix_t B = coco_block_matrix(dimension, block_sizes.data(), nb_blocks);
  coco_compute_blockrotation(B.rows(), rseed + 1000000, dimension, block_sizes.data(), nb_blocks);

  coco_compute_truncated_uniform_swap_permutation(P1.data(), rseed + 2000000, dimension, nb_swaps, swap_range);
  coco_compute_truncated_uniform_swap_permutation(P2.data(), rseed + 3000000, dimension, nb_swaps, swap_range);

  problem = f_ellipsoid_allocate(dimension, 1.0e6);
  problem = transform_vars_oscillate(std::move(problem));
  problem = transform_vars_permutation(std::move(problem), P2.data(), dimension);
  problem = transform_vars_blockrotation(std::move(problem), B.rows(), dimension, block_sizes.data(), nb_blocks);
  problem = transform_vars_permutation(std::move(problem), P1.data(), dimension);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);

  problem = transform_obj_norm_by_dim(std::move(problem));
  problem = transform_obj_shift(std::move(problem), fopt);

  coco_problem_set_id(problem.get(), problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem.get(), problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem.get(), "3-ill-conditioned");

  return problem;
}

/**
 * @brief Creates the ellipsoid problem for the constrained BBOB suite
 */
static coco_problem_ptr_t f_ellipsoid_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                 const size_t instance, const long rseed,
                                                                 const char *problem_id_template,
                                                                 const char *problem_name_template) {
  double fopt;
  coco_problem_ptr_t problem;

  std::vector<double> xopt(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt.data(), rseed, dimension);
  } else {
    bbob2009_compute_xopt(xopt.data(), rseed, dimension);
  }

  problem = f_ellipsoid_allocate(dimension, 1.0e6);

  /* TODO (NH): fopt -= problem->evaluate(all_zeros(dimension)) */
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);
  problem = transform_obj_shift(std::move(problem), fopt);

  coco_problem_set_id(problem.get(), problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem.get(), problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem.get(), "1-separable");

  return problem;
}

/**
 * @brief Creates the rotated ellipsoid problem for the constrained
 *        BBOB suite
 */
static coco_problem_ptr_t f_ellipsoid_rotated_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                         const size_t instance, const long rseed,
                                                                         const char *problem_id_template,
                                                                         const char *problem_name_template) {
  double fopt;
  coco_problem_ptr_t problem;

  std::vector<double> M(dimension * dimension);
  std::vector<double> b(dimension);

  std::vector<double> xopt(dimension);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt.data(), rseed, dimension);
  } else {
    bbob2009_compute_xopt(xopt.data(), rseed, dimension);
  }
  fopt = bbob2009_compute_fopt(function, instance);

  coco_matrix_t rot1(dimension, dimension);
  bbob2009_compute_rotation(rot1.rows(), rseed + 1000000, dimension);
  bbob2009_copy_rotation_matrix(rot1.rows(), M.data(), b.data(), dimension);

  problem = f_ellipsoid_allocate(dimension, 1.0e6);
  problem = transform_vars_affine(std::move(problem), M.data(), b.data(), dimension);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);
  problem = transform_obj_shift(std::move(problem), fopt);

  coco_problem_set_id(problem.get(), problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem.get(), problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem.get(), "3-ill-conditioned");

  return problem;
}
