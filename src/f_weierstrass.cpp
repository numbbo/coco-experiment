/**
 * @file f_weierstrass.cpp
 * @brief Implementation of the Weierstrass function and problem.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

/** @brief Number of summands in the Weierstrass problem. */
static constexpr size_t F_WEIERSTRASS_SUMMANDS = 12;

/**
 * @brief The Weierstrass function.
 */
class f_weierstrass_t final : public coco_problem_s {
public:
  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double raw(std::span<const double> x) const;

  explicit f_weierstrass_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {

    init_from_scalars("Weierstrass function", -5.0, 5.0, 0.);
    set_id("%s_d%02zu", "weierstrass", number_of_variables);

    f0 = 0.0;
    for (size_t i = 0; i < F_WEIERSTRASS_SUMMANDS; ++i) {
      ak[i] = pow(0.5, static_cast<double>(i));
      bk[i] = pow(3., static_cast<double>(i));
      f0 += ak[i] * cos(2 * coco_pi * bk[i] * 0.5);
    }

    /* Compute best solution */
    std::vector<double> non_unique_best_value(number_of_variables, 0.0);
    evaluate_function(non_unique_best_value, best_value);
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
  double f0{};
  double ak[F_WEIERSTRASS_SUMMANDS]{};
  double bk[F_WEIERSTRASS_SUMMANDS]{};
};

/**
 * @brief Implements the Weierstrass function without connections to any COCO
 * structures.
 */
double f_weierstrass_t::raw(std::span<const double> x) const {

  double result = 0.0;
  for (size_t i = 0; i < number_of_variables; ++i) {
    for (size_t j = 0; j < F_WEIERSTRASS_SUMMANDS; ++j) {
      result += cos(2 * coco_pi * (x[i] + 0.5) * bk[j]) * ak[j];
    }
  }
  result = 10.0 * pow(result / static_cast<double>(number_of_variables) - f0, 3.0);

  return result;
}

/**
 * @brief The BBOB Weierstrass function f16, also in the bbob-boxed suite.
 *
 * f(x) = f_weierstrass(R Λ^(1/100) Q T_osz(R (x - xopt))) + fopt + 10 / n f_pen(x).
 */
class bbob_weierstrass_t final : public bbob_problem_t<bbob_weierstrass_t> {
public:
  bbob_weierstrass_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                     const coco_bbob_variant_t variant, const char *problem_id_template,
                     const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "4-multi-modal", problem_id_template, problem_name_template),
        shift(bbob_xopt(rseed, dimension, variant)),
        rotation(bbob2009_compute_rotation(rseed + 1000000, dimension), std::vector<double>(dimension)),
        conditioned_rotation(conditioned_rotation_matrix(rseed, dimension), std::vector<double>(dimension)),
        weierstrass(dimension), fopt(bbob2009_compute_fopt(function, instance)),
        penalty_factor(10.0 / static_cast<double>(dimension)), z1(dimension), z2(dimension) {
    /* The optimum of the Weierstrass function is the origin, map it back through the transformations */
    conditioned_rotation.inverse(best_parameter, z1);
    oscillate.inverse(z1, z1);
    rotation.inverse(z1, z2);
    shift.inverse(z2, best_parameter);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = coco_boundary_penalty(x, smallest_values_of_interest, largest_values_of_interest);
    shift.apply(x, z1);
    rotation.apply(z1, z2);
    oscillate.apply(z2, z2);
    conditioned_rotation.apply(z2, z1);
    const double y = weierstrass.raw(z1) + fopt;
    return y + penalty_factor * penalty;
  }

private:
  /** @brief Returns R Λ^(1/100) Q with the rotations R and Q. */
  static coco_matrix_t conditioned_rotation_matrix(const long rseed, const size_t dimension) {
    const double condition = 100.0;
    const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
    const coco_matrix_t rot2 = bbob2009_compute_rotation(rseed, dimension);
    std::vector<double> scales(dimension);
    for (size_t k = 0; k < dimension; ++k) {
      const double base = 1.0 / sqrt(condition);
      const double exponent = 1.0 * static_cast<int>(k) / (static_cast<double>(dimension) - 1.0);
      scales[k] = pow(base, exponent);
    }
    return coco_matrix_product_with_diagonal(rot1, scales, rot2);
  }

  const vars_shift_t shift;
  const vars_affine_t rotation;
  const vars_oscillate_t oscillate{};
  const vars_affine_t conditioned_rotation;
  const f_weierstrass_t weierstrass;
  const double fopt;
  const double penalty_factor;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};

/**
 * @brief The Weierstrass function f16 of the bbob-largescale suite, with permuted block rotations.
 *
 * f(x) = f_weierstrass(R Λ^(1/100) Q T_osz(R (x - xopt))) + 10 / n f_pen(x) + fopt with the permuted block rotations
 * R and Q.
 */
class bbob_weierstrass_permblockdiag_t final : public bbob_problem_t<bbob_weierstrass_permblockdiag_t> {
public:
  bbob_weierstrass_permblockdiag_t(const size_t function, const size_t dimension, const size_t instance,
                                   const long rseed, const coco_bbob_variant_t variant, const char *problem_id_template,
                                   const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "4-multi-modal", problem_id_template, problem_name_template),
        shift(bbob_xopt(rseed, dimension, variant)),
        rotation1(dimension, rseed + 5000000, rseed + 1000000, rseed + 6000000),
        rotation2(dimension, rseed + 3000000, rseed, rseed + 4000000), conditioning(1.0 / 100.0, dimension),
        weierstrass(dimension), fopt(bbob2009_compute_fopt(function, instance)),
        penalty_factor(10.0 / static_cast<double>(dimension)), z1(dimension), z2(dimension) {
    /* The optimum of the Weierstrass function is the origin, map it back through the transformations */
    rotation1.inverse(best_parameter, z1);
    conditioning.inverse(z1, z1);
    rotation2.inverse(z1, z2);
    oscillate.inverse(z2, z2);
    rotation1.inverse(z2, z1);
    shift.inverse(z1, best_parameter);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = coco_boundary_penalty(x, smallest_values_of_interest, largest_values_of_interest);
    shift.apply(x, z1);
    rotation1.apply(z1, z2);
    oscillate.apply(z2, z2);
    rotation2.apply(z2, z1);
    conditioning.apply(z1, z1);
    rotation1.apply(z1, z2);
    const double y = weierstrass.raw(z2) + penalty_factor * penalty;
    return y + fopt;
  }

private:
  const vars_shift_t shift;
  const vars_permblockdiag_t rotation1;
  const vars_oscillate_t oscillate{};
  const vars_permblockdiag_t rotation2;
  const vars_conditioning_t conditioning;
  const f_weierstrass_t weierstrass;
  const double fopt;
  const double penalty_factor;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};
