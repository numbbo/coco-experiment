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

/**
 * @brief The ellipsoid function.
 */
struct f_ellipsoid_t final : public coco_problem_s {
  std::vector<double> scales; /**< @brief conditioning^(i / (n - 1)) */

  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double raw(const double *x) const;

  f_ellipsoid_t(const size_t number_of_variables, const double conditioning)
      : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("ellipsoid function", -5.0, 5.0, 0.0);
    set_id("%s_d%02zu", "ellipsoid", number_of_variables);

    scales.resize(number_of_variables);
    for (size_t i = 0; i < number_of_variables; ++i) {
      const double exponent = 1.0 * static_cast<double>(i) / (static_cast<double>(number_of_variables) - 1.0);
      scales[i] = pow(conditioning, exponent);
    }

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

  /**
   * @brief Evaluates the gradient of the ellipsoid function.
   */
  void evaluate_gradient(const double *x, double *y) override {
    for (size_t i = 0; i < number_of_variables; ++i)
      y[i] = 2.0 * scales[i] * x[i];
  }
};

/**
 * @brief Implements the ellipsoid function without connections to any COCO
 * structures.
 */
double f_ellipsoid_t::raw(const double *x) const {

  double result = x[0] * x[0];
  for (size_t i = 1; i < number_of_variables; ++i)
    result += scales[i] * x[i] * x[i];

  return result;
}

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
                                                            const coco_bbob_variant_t variant,
                                                            const char *problem_id_template,
                                                            const char *problem_name_template) {
  coco_problem_ptr_t problem;

  std::vector<double> xopt(dimension);
  double fopt = bbob2009_compute_fopt(function, instance);
  bbob_compute_xopt(xopt.data(), rseed, dimension, variant);

  problem = f_ellipsoid_allocate(dimension, 1.0e6);
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
 * @brief Creates the BBOB rotated ellipsoid problem.
 */
static coco_problem_ptr_t f_ellipsoid_rotated_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                    const size_t instance, const long rseed,
                                                                    const coco_bbob_variant_t variant,
                                                                    const double conditioning,
                                                                    const char *problem_id_template,
                                                                    const char *problem_name_template) {
  coco_problem_ptr_t problem;

  coco_matrix_t M(dimension, dimension);
  std::vector<double> b(dimension);

  std::vector<double> xopt(dimension);
  bbob_compute_xopt(xopt.data(), rseed, dimension, variant);
  double fopt = bbob2009_compute_fopt(function, instance);

  const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
  M = rot1;
  problem = f_ellipsoid_allocate(dimension, conditioning);
  problem = transform_vars_oscillate(std::move(problem));
  problem = transform_vars_affine(std::move(problem), M, b);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);
  problem = transform_obj_shift(std::move(problem), fopt);
  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("3-ill-conditioned");

  return problem;
}

/**
 * @brief Creates the BBOB permuted block-rotated ellipsoid problem.
 */
static coco_problem_ptr_t f_ellipsoid_permblockdiag_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                          const size_t instance, const long rseed,
                                                                          const coco_bbob_variant_t variant,
                                                                          const char *problem_id_template,
                                                                          const char *problem_name_template) {
  coco_problem_ptr_t problem;
  std::vector<size_t> P1(dimension);
  std::vector<size_t> P2(dimension);
  size_t nb_blocks;

  std::vector<size_t> block_sizes = coco_get_block_sizes(&nb_blocks, dimension, "bbob-largescale");
  size_t swap_range = coco_get_swap_range(dimension, "bbob-largescale");
  size_t nb_swaps = coco_get_nb_swaps(dimension, "bbob-largescale");

  std::vector<double> xopt(dimension);
  bbob_compute_xopt(xopt.data(), rseed, dimension, variant);
  double fopt = bbob2009_compute_fopt(function, instance);

  const coco_block_diagonal_matrix_t B = coco_compute_blockrotation(rseed + 1000000, block_sizes.data(), nb_blocks);

  coco_compute_truncated_uniform_swap_permutation(P1.data(), rseed + 2000000, dimension, nb_swaps, swap_range);
  coco_compute_truncated_uniform_swap_permutation(P2.data(), rseed + 3000000, dimension, nb_swaps, swap_range);

  problem = f_ellipsoid_allocate(dimension, 1.0e6);
  problem = transform_vars_oscillate(std::move(problem));
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

/**
 * @brief Creates the ellipsoid problem for the constrained BBOB suite
 */
static coco_problem_ptr_t f_ellipsoid_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                 const size_t instance, const long rseed,
                                                                 const coco_bbob_variant_t variant,
                                                                 const char *problem_id_template,
                                                                 const char *problem_name_template) {
  coco_problem_ptr_t problem;

  std::vector<double> xopt(dimension);
  double fopt = bbob2009_compute_fopt(function, instance);
  bbob_compute_xopt(xopt.data(), rseed, dimension, variant);

  problem = f_ellipsoid_allocate(dimension, 1.0e6);

  /* TODO (NH): fopt -= problem->evaluate(all_zeros(dimension)) */
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("1-separable");

  return problem;
}

/**
 * @brief Creates the rotated ellipsoid problem for the constrained
 *        BBOB suite
 */
static coco_problem_ptr_t f_ellipsoid_rotated_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                         const size_t instance, const long rseed,
                                                                         const coco_bbob_variant_t variant,
                                                                         const char *problem_id_template,
                                                                         const char *problem_name_template) {
  coco_problem_ptr_t problem;

  coco_matrix_t M(dimension, dimension);
  std::vector<double> b(dimension);

  std::vector<double> xopt(dimension);
  bbob_compute_xopt(xopt.data(), rseed, dimension, variant);
  double fopt = bbob2009_compute_fopt(function, instance);

  const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
  M = rot1;

  problem = f_ellipsoid_allocate(dimension, 1.0e6);
  problem = transform_vars_affine(std::move(problem), M, b);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("3-ill-conditioned");

  return problem;
}
