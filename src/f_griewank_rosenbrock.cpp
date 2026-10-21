/**
 * @file f_griewank_rosenbrock.cpp
 * @brief Implementation of the Griewank-Rosenbrock function and problem.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

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
 * @brief The BBOB Griewank-Rosenbrock function f19, also in the bbob-boxed and bbob-noisy suites.
 *
 * f(x) = f_griewank_rosenbrock(s R x + 0.5) + fopt with s = max(1, sqrt(n) / 8).
 */
class bbob_griewank_rosenbrock_t final : public bbob_problem_t<bbob_griewank_rosenbrock_t> {
public:
  bbob_griewank_rosenbrock_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                             const double facftrue, const char *problem_id_template, const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "4-multi-modal", problem_id_template, problem_name_template),
        rotation(scaled_rotation_matrix(rseed, dimension), std::vector<double>(dimension)),
        shift(std::vector<double>(dimension, -0.5)), griewank_rosenbrock(dimension, facftrue),
        fopt(bbob2009_compute_fopt(function, instance)), z1(dimension), z2(dimension) {
    /* The optimum of the Griewank-Rosenbrock function is (1, ..., 1), map it back through the transformations */
    const std::vector<double> ones(dimension, 1.0);
    shift.inverse(ones, z1);
    rotation.inverse(z1, best_parameter);
    /* Evaluating at the best parameter could add rounding errors */
    best_value[0] = griewank_rosenbrock.raw(ones) + fopt;
  }

  double operator()(std::span<const double> x) {
    rotation.apply(x, z1);
    shift.apply(z1, z2);
    return griewank_rosenbrock.raw(z2) + fopt;
  }

private:
  /** @brief Returns s R with the rotation R. */
  static coco_matrix_t scaled_rotation_matrix(const long rseed, const size_t dimension) {
    coco_matrix_t M = bbob2009_compute_rotation(rseed, dimension);
    const double scales = coco_double_max(1., sqrt(static_cast<double>(dimension)) / 8.);
    for (size_t i = 0; i < dimension; ++i) {
      for (size_t j = 0; j < dimension; ++j)
        M(i, j) *= scales;
    }
    return M;
  }

  const vars_affine_t rotation;
  const vars_shift_t shift; /**< @brief Adds 0.5 */
  const f_griewank_rosenbrock_t griewank_rosenbrock;
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};

/**
 * @brief The Griewank-Rosenbrock function f19 of the bbob-largescale suite, with a permuted block rotation.
 *
 * f(x) = f_griewank_rosenbrock(s R x + 0.5) + fopt with the permuted block rotation R and
 * s = max(1, sqrt(block size) / 8).
 */
class bbob_griewank_rosenbrock_permblockdiag_t final : public bbob_problem_t<bbob_griewank_rosenbrock_permblockdiag_t> {
public:
  bbob_griewank_rosenbrock_permblockdiag_t(const size_t function, const size_t dimension, const size_t instance,
                                           const long rseed, const double facftrue, const char *problem_id_template,
                                           const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "4-multi-modal", problem_id_template, problem_name_template),
        rotation(dimension, rseed + 2000000, rseed, rseed + 3000000),
        scale(coco_double_max(1.0, sqrt(static_cast<double>(coco_rotation_matrix_block_size(dimension))) / 8.0),
              dimension),
        shift(std::vector<double>(dimension, -0.5)), griewank_rosenbrock(dimension, facftrue),
        fopt(bbob2009_compute_fopt(function, instance)), z1(dimension), z2(dimension) {
    /* The optimum of the Griewank-Rosenbrock function is (1, ..., 1), map it back through the transformations */
    const std::vector<double> ones(dimension, 1.0);
    shift.inverse(ones, z1);
    scale.inverse(z1, z1);
    rotation.inverse(z1, best_parameter);
    /* Evaluating at the best parameter could add rounding errors */
    best_value[0] = griewank_rosenbrock.raw(ones) + fopt;
  }

  double operator()(std::span<const double> x) {
    rotation.apply(x, z1);
    scale.apply(z1, z1);
    shift.apply(z1, z2);
    return griewank_rosenbrock.raw(z2) + fopt;
  }

private:
  const vars_permblockdiag_t rotation;
  const vars_scale_t scale;
  const vars_shift_t shift; /**< @brief Adds 0.5 */
  const f_griewank_rosenbrock_t griewank_rosenbrock;
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};
