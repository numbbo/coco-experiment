/**
 * @file f_schaffers.cpp
 * @brief Implementation of the Schaffer's F7 function and problem,
 * transformations not implemented for the moment.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

/**
 * @brief Implements the Schaffer's F7 function without connections to any COCO
 * structures.
 */
static double f_schaffers_raw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  size_t i = 0;

  assert(number_of_variables > 1);

  /* Computation core */
  double result = 0.0;
  for (i = 0; i < number_of_variables - 1; ++i) {
    const double tmp = x[i] * x[i] + x[i + 1] * x[i + 1];
    if (coco_is_inf(tmp) && coco_is_nan(sin(50.0 * pow(tmp, 0.1)))) /* sin(inf) -> nan */
      /* the second condition is necessary to pass the integration tests under Windows and Linux */
      return tmp;
    result += pow(tmp, 0.25) * (1.0 + pow(sin(50.0 * pow(tmp, 0.1)), 2.0));
  }
  result = pow(result / (static_cast<double>(number_of_variables) - 1.0), 2.0);

  return result;
}

/**
 * @brief The BBOB Schaffer's F7 functions f17 and f18, also in the bbob-boxed and bbob-noisy suites.
 *
 * f(x) = f_schaffers(Λ^conditioning Q T_asy^0.5(R (x - xopt))) + fopt + penalty_scale f_pen(x).
 */
class bbob_schaffers_t final : public bbob_problem_t<bbob_schaffers_t> {
public:
  bbob_schaffers_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                   const coco_bbob_variant_t variant, const double conditioning, const double penalty_scale,
                   const char *problem_id_template, const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "4-multi-modal", problem_id_template, problem_name_template),
        shift(bbob_xopt(rseed, dimension, variant)),
        rotation(bbob2009_compute_rotation(rseed + 1000000, dimension), std::vector<double>(dimension)),
        asymmetric(0.5, dimension), conditioned_rotation(conditioned_rotation_matrix(rseed, dimension, conditioning),
                                                         std::vector<double>(dimension)),
        fopt(bbob2009_compute_fopt(function, instance)), penalty_scale(penalty_scale), z1(dimension), z2(dimension) {
    /* The optimum of Schaffer's function is the origin, map it back through the transformations */
    conditioned_rotation.inverse(best_parameter, z1);
    asymmetric.inverse(z1, z1);
    rotation.inverse(z1, z2);
    shift.inverse(z2, best_parameter);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = coco_boundary_penalty(x, smallest_values_of_interest, largest_values_of_interest);
    shift.apply(x, z1);
    rotation.apply(z1, z2);
    asymmetric.apply(z2, z2);
    conditioned_rotation.apply(z2, z1);
    const double y = f_schaffers_raw(z1) + fopt;
    return y + penalty_scale * penalty;
  }

private:
  /** @brief Returns Λ^conditioning Q with the rotation Q. */
  static coco_matrix_t conditioned_rotation_matrix(const long rseed, const size_t dimension,
                                                   const double conditioning) {
    coco_matrix_t M = bbob2009_compute_rotation(rseed, dimension);
    for (size_t i = 0; i < dimension; ++i) {
      for (size_t j = 0; j < dimension; ++j) {
        double exponent = 1.0 * static_cast<int>(i) / (static_cast<double>(dimension) - 1.0);
        M(i, j) = M(i, j) * pow(sqrt(conditioning), exponent);
      }
    }
    return M;
  }

  const vars_shift_t shift;
  const vars_affine_t rotation;
  const vars_asymmetric_t asymmetric;
  const vars_affine_t conditioned_rotation;
  const double fopt;
  const double penalty_scale;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};

/**
 * @brief The Schaffer's F7 functions f17 and f18 of the bbob-largescale suite, with permuted block rotations.
 *
 * f(x) = f_schaffers(Λ^conditioning Q T_asy^0.5(R (x - xopt))) + 10 f_pen(x) + fopt with the permuted block rotations
 * R and Q.
 */
class bbob_schaffers_permblockdiag_t final : public bbob_problem_t<bbob_schaffers_permblockdiag_t> {
public:
  bbob_schaffers_permblockdiag_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                                 const coco_bbob_variant_t variant, const double condition,
                                 const char *problem_id_template, const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "4-multi-modal", problem_id_template, problem_name_template),
        shift(bbob_xopt(rseed, dimension, variant)),
        rotation1(dimension, rseed + 5000000, rseed + 1000000, rseed + 6000000), asymmetric(0.5, dimension),
        rotation2(dimension, rseed + 3000000, rseed, rseed + 4000000), conditioning(condition, dimension),
        fopt(bbob2009_compute_fopt(function, instance)), z1(dimension), z2(dimension) {
    /* The optimum of Schaffer's function is the origin, map it back through the transformations */
    conditioning.inverse(best_parameter, z1);
    rotation2.inverse(z1, z2);
    asymmetric.inverse(z2, z2);
    rotation1.inverse(z2, z1);
    shift.inverse(z1, best_parameter);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = coco_boundary_penalty(x, smallest_values_of_interest, largest_values_of_interest);
    shift.apply(x, z1);
    rotation1.apply(z1, z2);
    asymmetric.apply(z2, z2);
    rotation2.apply(z2, z1);
    conditioning.apply(z1, z1);
    const double y = f_schaffers_raw(z1) + penalty_factor * penalty;
    return y + fopt;
  }

private:
  static constexpr double penalty_factor = 10.0;

  const vars_shift_t shift;
  const vars_permblockdiag_t rotation1;
  const vars_asymmetric_t asymmetric;
  const vars_permblockdiag_t rotation2;
  const vars_conditioning_t conditioning;
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};
