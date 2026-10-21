/**
 * @file f_rosenbrock.cpp
 * @brief Implementation of the Rosenbrock function and problem.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

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
 * @brief The location of the optimum of the BBOB Rosenbrock functions before the shift to (1, ..., 1): 0.75 xopt.
 */
static std::vector<double> rosenbrock_xopt(const long rseed, const size_t dimension,
                                           const coco_bbob_variant_t variant) {
  std::vector<double> xopt = bbob_xopt(rseed, dimension, variant);
  for (size_t i = 0; i < dimension; ++i)
    xopt[i] *= 0.75;
  return xopt;
}

/**
 * @brief The scaling factor of the BBOB Rosenbrock functions: max(1, sqrt(n) / 8), where n is the block size in
 * the large-scale suite.
 */
static double rosenbrock_factor(const size_t dimension, const coco_bbob_variant_t variant) {
  if (variant == coco_bbob_variant_t::largescale) {
    const size_t block_size = coco_rotation_matrix_block_size(dimension);
    return coco_double_max(1.0, sqrt(static_cast<double>(block_size)) / 8.0);
  }
  return coco_double_max(1.0, sqrt(static_cast<double>(dimension)) / 8.0);
}

/**
 * @brief The BBOB Rosenbrock function f8, also in the bbob-boxed, bbob-largescale and bbob-noisy suites.
 *
 * f(x) = f_rosenbrock(s (x - 0.75 xopt) + 1) + fopt with s = max(1, sqrt(n) / 8), where n is the block size in the
 * large-scale suite; normalized by the dimension before fopt is added in the large-scale suite.
 */
class bbob_rosenbrock_t final : public bbob_problem_t<bbob_rosenbrock_t> {
public:
  bbob_rosenbrock_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                    const coco_bbob_variant_t variant, const char *problem_id_template,
                    const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "2-moderate", problem_id_template, problem_name_template),
        shift(rosenbrock_xopt(rseed, dimension, variant)), scale(rosenbrock_factor(dimension, variant), dimension),
        shift_to_one(std::vector<double>(dimension, -1.0)), fopt(bbob2009_compute_fopt(function, instance)),
        normalize(variant == coco_bbob_variant_t::largescale), z(dimension) {
    /* The optimum of the Rosenbrock function is (1, ..., 1), map it back through the transformations */
    std::fill(z.begin(), z.end(), 1.0);
    shift_to_one.inverse(z, z);
    scale.inverse(z, z);
    shift.inverse(z, best_parameter);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    shift.apply(x, z);
    scale.apply(z, z);
    shift_to_one.apply(z, z);
    double y = f_rosenbrock_raw(z);
    if (normalize)
      y *= coco_norm_by_dim_factor(number_of_variables);
    return y + fopt;
  }

private:
  const vars_shift_t shift;
  const vars_scale_t scale;
  const vars_shift_t shift_to_one; /**< @brief Adds 1 */
  const double fopt;
  const bool normalize;
  std::vector<double> z; /**< @brief The transformed decision values */
};

/**
 * @brief The BBOB rotated Rosenbrock function f9, also in the bbob-boxed suite.
 *
 * f(x) = f_rosenbrock(s R x + 0.5) + fopt with s = max(1, sqrt(n) / 8).
 */
class bbob_rosenbrock_rotated_t final : public bbob_problem_t<bbob_rosenbrock_rotated_t> {
public:
  bbob_rosenbrock_rotated_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                            const char *problem_id_template, const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "2-moderate", problem_id_template, problem_name_template),
        affine(scaled_rotation_matrix(rseed, dimension), std::vector<double>(dimension, 0.5)),
        fopt(bbob2009_compute_fopt(function, instance)), z(dimension) {
    /* The optimum of the Rosenbrock function is (1, ..., 1), map it back through the transformation */
    const std::vector<double> ones(dimension, 1.0);
    affine.inverse(ones, best_parameter);
    /* Evaluating at the best parameter could add rounding errors */
    best_value[0] = f_rosenbrock_raw(ones) + fopt;
  }

  double operator()(std::span<const double> x) {
    affine.apply(x, z);
    return f_rosenbrock_raw(z) + fopt;
  }

private:
  /** @brief Returns s R with the rotation R. */
  static coco_matrix_t scaled_rotation_matrix(const long rseed, const size_t dimension) {
    const double factor = coco_double_max(1.0, sqrt(static_cast<double>(dimension)) / 8.0);
    coco_matrix_t M = bbob2009_compute_rotation(rseed, dimension);
    for (size_t row = 0; row < dimension; ++row) {
      for (size_t column = 0; column < dimension; ++column)
        M(row, column) = factor * M(row, column);
    }
    return M;
  }

  const vars_affine_t affine;
  const double fopt;
  std::vector<double> z; /**< @brief The transformed decision values */
};

/**
 * @brief The rotated Rosenbrock function f9 of the bbob-largescale suite, with a permuted block rotation.
 *
 * f(x) = f_rosenbrock(s R (x - 0.75 xopt) + 1) / n + fopt with the permuted block rotation R and
 * s = max(1, sqrt(block size) / 8).
 */
class bbob_rosenbrock_permblockdiag_t final : public bbob_problem_t<bbob_rosenbrock_permblockdiag_t> {
public:
  bbob_rosenbrock_permblockdiag_t(const size_t function, const size_t dimension, const size_t instance,
                                  const long rseed, const coco_bbob_variant_t variant, const char *problem_id_template,
                                  const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "2-moderate", problem_id_template, problem_name_template),
        shift(rosenbrock_xopt(rseed, dimension, variant)), rotation(dimension, rseed + 2000000, rseed, rseed + 3000000),
        scale(rosenbrock_factor(dimension, coco_bbob_variant_t::largescale), dimension),
        shift_to_one(std::vector<double>(dimension, -1.0)), fopt(bbob2009_compute_fopt(function, instance)),
        z1(dimension), z2(dimension) {
    /* The optimum of the Rosenbrock function is (1, ..., 1), map it back through the transformations */
    std::fill(z1.begin(), z1.end(), 1.0);
    shift_to_one.inverse(z1, z1);
    scale.inverse(z1, z1);
    rotation.inverse(z1, z2);
    shift.inverse(z2, best_parameter);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    shift.apply(x, z1);
    rotation.apply(z1, z2);
    scale.apply(z2, z2);
    shift_to_one.apply(z2, z2);
    return f_rosenbrock_raw(z2) * coco_norm_by_dim_factor(number_of_variables) + fopt;
  }

private:
  const vars_shift_t shift;
  const vars_permblockdiag_t rotation;
  const vars_scale_t scale;
  const vars_shift_t shift_to_one; /**< @brief Adds 1 */
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};
