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
 * @brief The Griewank-Rosenbrock function.
 */
class f_griewank_rosenbrock_t final : public coco_problem_s {
public:
  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double raw(std::span<const double> x) const;

  f_griewank_rosenbrock_t(const size_t number_of_variables, double facftrue)
      : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("Griewank Rosenbrock function", -5.0, 5.0, 1);
    set_id("%s_d%02zu", "griewank_rosenbrock", number_of_variables);

    this->facftrue = facftrue;

    /* Compute best solution */
    evaluate_function(best_parameter, best_value);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives == 1);
    y[0] = raw(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

private:
  double facftrue{};
};

/**
 * @brief Implements the Griewank-Rosenbrock function without connections to any
 * COCO structures.
 */
double f_griewank_rosenbrock_t::raw(std::span<const double> x) const {

  size_t i = 0;
  double tmp = 0;

  /* Computation core */
  double result = 0.0;
  for (i = 0; i < number_of_variables - 1; ++i) {
    const double c1 = x[i] * x[i] - x[i + 1];
    const double c2 = 1.0 - x[i];
    tmp = 100.0 * c1 * c1 + c2 * c2;
    result += tmp / 4000. - cos(tmp);
  }
  result = facftrue + facftrue * result / static_cast<double>(number_of_variables - 1);

  return result;
}

/**
 * @brief Allocates the basic Griewank-Rosenbrock problem.
 */
static coco_problem_ptr_t f_griewank_rosenbrock_allocate(const size_t number_of_variables, double facftrue) {
  return std::make_unique<f_griewank_rosenbrock_t>(number_of_variables, facftrue);
}

/**
 * @brief Creates the BBOB Griewank-Rosenbrock problem.
 */
static coco_problem_ptr_t f_griewank_rosenbrock_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                      const size_t instance, const long rseed,
                                                                      const double facftrue,
                                                                      const char *problem_id_template,
                                                                      const char *problem_name_template) {
  coco_problem_ptr_t problem;

  coco_matrix_t M(dimension, dimension);
  std::vector<double> b(dimension);
  std::vector<double> shift(dimension);

  double fopt = bbob2009_compute_fopt(function, instance);
  for (size_t i = 0; i < dimension; ++i) {
    shift[i] = -0.5;
  }

  coco_matrix_t rot1 = bbob2009_compute_rotation(rseed, dimension);
  double scales = coco_double_max(1., sqrt(static_cast<double>(dimension)) / 8.);
  for (size_t i = 0; i < dimension; ++i) {
    for (size_t j = 0; j < dimension; ++j) {
      rot1(i, j) *= scales;
    }
  }

  problem = f_griewank_rosenbrock_allocate(dimension, facftrue);
  problem = transform_obj_shift(std::move(problem), fopt);
  problem = transform_vars_shift(std::move(problem), shift, 0);
  M = rot1;

  problem = transform_vars_affine(std::move(problem), M, b);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("4-multi-modal");

  return problem;
}

/**
 * @brief Creates the BBOB permuted block-rotated Griewank-Rosenbrock problem.
 */
static coco_problem_ptr_t f_griewank_rosenbrock_permblockdiag_bbob_bbob_problem_allocate(
    const size_t function, const size_t dimension, const size_t instance, const long rseed, const double facftrue,
    const char *problem_id_template, const char *problem_name_template) {
  coco_problem_ptr_t problem;

  std::vector<size_t> P1(dimension);
  std::vector<size_t> P2(dimension);
  size_t nb_blocks;

  std::vector<size_t> block_sizes = coco_get_block_sizes(&nb_blocks, dimension, "bbob-largescale");
  size_t block_size = coco_rotation_matrix_block_size(dimension);
  size_t swap_range = coco_get_swap_range(dimension, "bbob-largescale");
  size_t nb_swaps = coco_get_nb_swaps(dimension, "bbob-largescale");

  double fopt = bbob2009_compute_fopt(function, instance);
  double scales = coco_double_max(1.0, sqrt(static_cast<double>(block_size)) / 8.0);
  std::vector<double> shift(dimension);
  for (size_t i = 0; i < dimension; ++i) {
    shift[i] = -0.5;
  }

  const coco_block_diagonal_matrix_t B = coco_compute_blockrotation(rseed, block_sizes);
  coco_compute_truncated_uniform_swap_permutation(P1, rseed + 2000000, nb_swaps, swap_range);
  coco_compute_truncated_uniform_swap_permutation(P2, rseed + 3000000, nb_swaps, swap_range);

  problem = f_griewank_rosenbrock_allocate(dimension, facftrue);
  problem = transform_vars_shift(std::move(problem), shift, 0);
  problem = transform_vars_scale(std::move(problem), scales);
  problem = transform_vars_permutation(std::move(problem), P2);
  problem = transform_vars_blockrotation(std::move(problem), B);
  problem = transform_vars_permutation(std::move(problem), P1);

  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("4-multi-modal");

  return problem;
}
