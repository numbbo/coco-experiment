/**
 * @file f_katsuura.cpp
 * @brief Implementation of the Katsuura function and problem.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "coco_utilities.cpp"
#include "suite_bbob_legacy_code.cpp"

/**
 * @brief Implements the Katsuura function without connections to any COCO
 * structures.
 */
static double f_katsuura_raw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  const double exponent = 10. / pow(static_cast<double>(number_of_variables), 1.2);

  /* Computation core */
  double result = 1.0;
  for (size_t i = 0; i < number_of_variables; ++i) {
    double tmp = 0;
    for (int j = 1; j < 33; ++j) {
      const double tmp2 = ldexp(1.0, j); /* 2^j */
      tmp += fabs(tmp2 * x[i] - coco_double_round(tmp2 * x[i])) / tmp2;
    }
    tmp = 1.0 + (static_cast<double>(i) + 1) * tmp;
    result *= pow(tmp, exponent);
  }
  result =
      10. / (static_cast<double>(number_of_variables)) / (static_cast<double>(number_of_variables)) * (-1. + result);

  return result;
}

/**
 * @brief The BBOB Katsuura function f23, also in the bbob-boxed suite.
 *
 * f(x) = f_katsuura(Q Λ^100 R (x - xopt)) + fopt + f_pen(x).
 */
class bbob_katsuura_t final : public bbob_problem_t<bbob_katsuura_t> {
public:
  bbob_katsuura_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                  const coco_bbob_variant_t variant, const char *problem_id_template, const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "5-weakly-structured", problem_id_template,
                       problem_name_template),
        shift(bbob_xopt(rseed, dimension, variant)),
        conditioned_rotation(conditioned_rotation_matrix(rseed, dimension), std::vector<double>(dimension)),
        fopt(bbob2009_compute_fopt(function, instance)), z1(dimension), z2(dimension) {
    /* The optimum of the Katsuura function is the origin, map it back through the transformations */
    conditioned_rotation.inverse(best_parameter, z1);
    shift.inverse(z1, best_parameter);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = coco_boundary_penalty(x, smallest_values_of_interest, largest_values_of_interest);
    shift.apply(x, z1);
    conditioned_rotation.apply(z1, z2);
    /* There is no shift 'fopt' in the definition */
    const double y = f_katsuura_raw(z2) + fopt;
    return y + penalty_factor * penalty;
  }

private:
  /** @brief Returns Q Λ^100 R with the rotations Q and R. */
  static coco_matrix_t conditioned_rotation_matrix(const long rseed, const size_t dimension) {
    const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
    const coco_matrix_t rot2 = bbob2009_compute_rotation(rseed, dimension);
    std::vector<double> scales(dimension);
    for (size_t k = 0; k < dimension; ++k) {
      const double exponent = 1.0 * static_cast<int>(k) / (static_cast<double>(dimension) - 1.0);
      scales[k] = pow(sqrt(100), exponent);
    }
    return coco_matrix_product_with_diagonal(rot1, scales, rot2);
  }

  static constexpr double penalty_factor = 1.0;

  const vars_shift_t shift;
  const vars_affine_t conditioned_rotation;
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};

/**
 * @brief The Katsuura function f23 of the bbob-largescale suite, with permuted block rotations.
 *
 * f(x) = f_katsuura(Q Λ^100 R (x - xopt)) + f_pen(x) + fopt with the permuted block rotations R and Q.
 */
class bbob_katsuura_permblockdiag_t final : public bbob_problem_t<bbob_katsuura_permblockdiag_t> {
public:
  bbob_katsuura_permblockdiag_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                                const coco_bbob_variant_t variant, const char *problem_id_template,
                                const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "5-weakly-structured", problem_id_template,
                       problem_name_template),
        shift(bbob_xopt(rseed, dimension, variant)), rotation1(dimension, rseed + 3000000, rseed, rseed + 4000000),
        conditioning(100.0, dimension), rotation2(dimension, rseed + 5000000, rseed + 1000000, rseed + 6000000),
        fopt(bbob2009_compute_fopt(function, instance)), z1(dimension), z2(dimension) {
    /* The optimum of the Katsuura function is the origin, map it back through the transformations */
    rotation2.inverse(best_parameter, z1);
    conditioning.inverse(z1, z1);
    rotation1.inverse(z1, z2);
    shift.inverse(z2, best_parameter);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = coco_boundary_penalty(x, smallest_values_of_interest, largest_values_of_interest);
    shift.apply(x, z1);
    rotation1.apply(z1, z2);
    conditioning.apply(z2, z2);
    rotation2.apply(z2, z1);
    const double y = f_katsuura_raw(z1) + penalty_factor * penalty;
    /* TODO: documentation, there is no fopt in the definition of this function */
    return y + fopt;
  }

private:
  static constexpr double penalty_factor = 1.0;

  const vars_shift_t shift;
  const vars_permblockdiag_t rotation1;
  const vars_conditioning_t conditioning;
  const vars_permblockdiag_t rotation2;
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};
