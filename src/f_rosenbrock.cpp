/**
 * @file f_rosenbrock.cpp
 * @brief Implementation of the Rosenbrock function and problem.
 */

#include <assert.h>

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
 * @brief Implements the Rosenbrock function without connections to any COCO
 * structures.
 */
static double f_rosenbrock_raw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  size_t i = 0;
  double s1 = 0.0, s2 = 0.0, tmp;

  assert(number_of_variables > 1);

  for (i = 0; i < number_of_variables - 1; ++i) {
    tmp = (x[i] * x[i] - x[i + 1]);
    s1 += tmp * tmp;
    tmp = (x[i] - 1.0);
    s2 += tmp * tmp;
  }
  double result = 100.0 * s1 + s2;

  return result;
}

/**
 * @brief The Rosenbrock function.
 */
class f_rosenbrock_t final : public coco_problem_s {
public:
  explicit f_rosenbrock_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("Rosenbrock function", -5.0, 5.0, 1.0);
    set_id("%s_d%02zu", "rosenbrock", number_of_variables);

    /* Compute best solution */
    evaluate_function(best_parameter, best_value);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives == 1);
    y[0] = f_rosenbrock_raw(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Allocates the basic Rosenbrock problem.
 */
static coco_problem_ptr_t f_rosenbrock_allocate(const size_t number_of_variables) {
  return std::make_unique<f_rosenbrock_t>(number_of_variables);
}

/**
 * @brief Creates the BBOB Rosenbrock problem.
 */
static coco_problem_ptr_t f_rosenbrock_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                             const size_t instance, const long rseed,
                                                             const coco_bbob_variant_t variant,
                                                             const char *problem_id_template,
                                                             const char *problem_name_template) {
  coco_problem_ptr_t problem;
  size_t block_size;
  double factor;

  std::vector<double> minus_one(dimension);
  std::vector<double> xopt(dimension);
  bbob_compute_xopt(xopt, rseed, variant);
  for (size_t i = 0; i < dimension; ++i) {
    minus_one[i] = -1.0;
    xopt[i] *= 0.75;
  }
  double fopt = bbob2009_compute_fopt(function, instance);
  if (variant == coco_bbob_variant_t::largescale) {
    block_size = coco_rotation_matrix_block_size(dimension);
    factor = coco_double_max(1.0, sqrt(static_cast<double>(block_size)) / 8.0);
  } else {
    factor = coco_double_max(1.0, sqrt(static_cast<double>(dimension)) / 8.0);
  }

  problem = f_rosenbrock_allocate(dimension);
  problem = transform_vars_shift(std::move(problem), minus_one, 0);
  problem = transform_vars_scale(std::move(problem), factor);
  problem = transform_vars_shift(std::move(problem), xopt, 0);

  /*if large scale test-bed, normalize by dim*/
  if (variant == coco_bbob_variant_t::largescale) {
    problem = transform_obj_norm_by_dim(std::move(problem));
  }

  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("2-moderate");

  return problem;
}

/**
 * @brief Creates the BBOB rotated Rosenbrock problem.
 */
static coco_problem_ptr_t f_rosenbrock_rotated_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                     const size_t instance, const long rseed,
                                                                     const char *problem_id_template,
                                                                     const char *problem_name_template) {

  coco_problem_ptr_t problem;
  coco_matrix_t M(dimension, dimension);
  std::vector<double> b(dimension);

  double fopt = bbob2009_compute_fopt(function, instance);
  const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed, dimension);

  double factor = coco_double_max(1.0, sqrt(static_cast<double>(dimension)) / 8.0);
  /* Compute affine transformation */
  for (size_t row = 0; row < dimension; ++row) {
    for (size_t column = 0; column < dimension; ++column) {
      M(row, column) = factor * rot1(row, column);
    }
    b[row] = 0.5;
  }

  problem = f_rosenbrock_allocate(dimension);
  problem = transform_vars_affine(std::move(problem), M, b);
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("2-moderate");

  return problem;
}

/**
 * @brief Creates the BBOB permuted block-rotated Rosenbrock problem.
 */
static coco_problem_ptr_t
f_rosenbrock_permblockdiag_bbob_problem_allocate(const size_t function, const size_t dimension, const size_t instance,
                                                 const long rseed, const coco_bbob_variant_t variant,
                                                 const char *problem_id_template, const char *problem_name_template) {

  coco_problem_ptr_t problem;

  std::vector<size_t> P1(dimension);
  std::vector<size_t> P2(dimension);
  size_t nb_blocks;
  std::vector<double> best_parameter(dimension); /* Manh: will serve to set the optimal solution "manually"*/

  std::vector<size_t> block_sizes = coco_get_block_sizes(&nb_blocks, dimension, "bbob-largescale");
  size_t block_size = coco_rotation_matrix_block_size(dimension);

  size_t swap_range = coco_get_swap_range(dimension, "bbob-largescale");
  size_t nb_swaps = coco_get_nb_swaps(dimension, "bbob-largescale");

  double fopt = bbob2009_compute_fopt(function, instance);
  double factor = coco_double_max(1.0, sqrt(static_cast<double>(block_size)) / 8.0);
  std::vector<double> minus_one(dimension);
  std::vector<double> xopt(dimension);
  bbob_compute_xopt(xopt, rseed, variant);
  for (size_t i = 0; i < dimension; ++i) {
    minus_one[i] = -1.0;
    xopt[i] *= 0.75;
  }

  const coco_block_diagonal_matrix_t B = coco_compute_blockrotation(rseed, block_sizes);
  coco_compute_truncated_uniform_swap_permutation(P1, rseed + 2000000, nb_swaps, swap_range);
  coco_compute_truncated_uniform_swap_permutation(P2, rseed + 3000000, nb_swaps, swap_range);

  problem = f_rosenbrock_allocate(dimension);
  problem = transform_vars_shift(std::move(problem), minus_one, 0);
  problem = transform_vars_scale(std::move(problem), factor);
  problem = transform_vars_permutation(std::move(problem), P2);
  problem = transform_vars_blockrotation(std::move(problem), B);
  problem = transform_vars_permutation(std::move(problem), P1);

  problem = transform_vars_shift(std::move(problem), xopt, 0);
  problem = transform_obj_norm_by_dim(std::move(problem));
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("2-moderate");

  return problem;
}
