/**
 * @file f_step_ellipsoid.cpp
 * @brief Implementation of the step ellipsoid function and problem.
 *
 * The BBOB step ellipsoid function intertwines the variable and objective
 * transformations in such a way that it is hard to devise a composition of
 * generic transformations to implement it. In the end one would have to
 * implement several custom transformations which would be used solely by this
 * problem. Therefore we opt to implement it as a monolithic function instead.
 *
 * TODO: It would be nice to have a generic step ellipsoid function to
 * complement this one.
 */
#include <assert.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "coco_utilities.cpp"
#include "suite_bbob_legacy_code.cpp"

#include "bbob_problem.cpp"
#include "transform_functions.cpp"

/**
 * @brief The BBOB step ellipsoid function f7, also in the bbob-boxed and bbob-noisy suites.
 *
 * f(x) = 0.1 max(|ẑ_1| 1e-4, f_ellipsoid^100(Q round(ẑ))) + penalty_scale f_pen(x) + fopt with
 * ẑ = Λ^10 R (x - xopt) and the rounding to steps of 1 or 0.1 (for |ẑ_i| <= 0.5).
 */
class bbob_step_ellipsoid_t final : public bbob_problem_t<bbob_step_ellipsoid_t> {
public:
  bbob_step_ellipsoid_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                        const coco_bbob_variant_t variant, const double penalty_scale, const char *problem_id_template,
                        const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "2-moderate", problem_id_template, problem_name_template),
        shift(bbob_xopt(rseed, dimension, variant)),
        conditioned_rotation(conditioned_rotation_matrix(rseed, dimension), std::vector<double>(dimension)),
        round_step(10.0),
        rotation(bbob2009_compute_rotation(rseed + 1000000, dimension), std::vector<double>(dimension)),
        scales(dimension), fopt(bbob2009_compute_fopt(function, instance)), penalty_scale(penalty_scale), z1(dimension),
        z2(dimension) {
    static const double condition = 100;
    assert(dimension > 1);
    for (size_t i = 0; i < dimension; ++i)
      scales[i] = pow(condition, static_cast<double>(i) / (static_cast<double>(dimension) - 1.0));

    /* The optimum is at xopt, where ẑ = 0 */
    shift.inverse(best_parameter, best_parameter);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty =
        penalty_scale * coco_boundary_penalty(x, smallest_values_of_interest, largest_values_of_interest);
    shift.apply(x, z1);
    conditioned_rotation.apply(z1, z2);
    const double zhat_1 = z2[0];
    round_step.apply(z2, z2);
    rotation.apply(z2, z1);

    double y = 0.0;
    for (size_t i = 0; i < number_of_variables; ++i)
      y += scales[i] * z1[i] * z1[i];
    return 0.1 * coco_double_max(fabs(zhat_1) * 1.0e-4, y) + penalty + fopt;
  }

private:
  /** @brief Returns Λ^10 R with the rotation R. */
  static coco_matrix_t conditioned_rotation_matrix(const long rseed, const size_t dimension) {
    static const double condition = 100;
    coco_matrix_t M = bbob2009_compute_rotation(rseed, dimension);
    for (size_t i = 0; i < dimension; ++i) {
      const double c1 = sqrt(pow(condition / 10., static_cast<double>(i) / static_cast<double>(dimension - 1)));
      for (size_t j = 0; j < dimension; ++j)
        M(i, j) = c1 * M(i, j);
    }
    return M;
  }

  const vars_shift_t shift;
  const vars_affine_t conditioned_rotation;
  const vars_round_step_t round_step;
  const vars_affine_t rotation;
  std::vector<double> scales; /**< @brief condition^(i / (n - 1)) */
  const double fopt;
  const double penalty_scale;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};

/**
 * @brief The step ellipsoid function f7 of the bbob-largescale suite, with permuted block rotations.
 *
 * f(x) = 0.1 max(|ẑ_1| n 1e-4, f_ellipsoid^100(Q round(ẑ))) / n + f_pen(x) + fopt with ẑ = Λ^10 R (x - xopt), the
 * rounding to steps of 1 or 0.1 (for |ẑ_i| <= 0.5), and the permuted block rotations R and Q.
 */
class bbob_step_ellipsoid_permblockdiag_t final : public bbob_problem_t<bbob_step_ellipsoid_permblockdiag_t> {
public:
  bbob_step_ellipsoid_permblockdiag_t(const size_t function, const size_t dimension, const size_t instance,
                                      const long rseed, const coco_bbob_variant_t variant,
                                      const char *problem_id_template, const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "2-moderate", problem_id_template, problem_name_template),
        shift(bbob_xopt(rseed, dimension, variant)), rotation1(dimension, rseed + 2000000, rseed, rseed + 3000000),
        conditioning(10.0, dimension), round_step(10.0),
        rotation2(dimension, rseed + 4000000, rseed + 1000000, rseed + 5000000), scales(dimension),
        fopt(bbob2009_compute_fopt(function, instance)), z1(dimension), z2(dimension) {
    static const double condition = 100;
    for (size_t i = 0; i < dimension; ++i)
      scales[i] = pow(condition, static_cast<double>(i) / (static_cast<double>(dimension) - 1.0));

    /* The optimum of the ellipsoid is the origin, map it back through the transformations; the rounding keeps it */
    rotation2.inverse(best_parameter, z1);
    round_step.apply(z1, z1);
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
    /* multiplied by the dimension to counter-balance the normalization */
    const double zhat_1 = fabs(z2[0]) * static_cast<double>(number_of_variables);
    round_step.apply(z2, z2);
    rotation2.apply(z2, z1);

    double y = 0.0;
    for (size_t i = 0; i < number_of_variables; ++i)
      y += scales[i] * z1[i] * z1[i];
    y = 0.1 * coco_double_max(zhat_1 * 1.0e-4, y);
    return y * coco_norm_by_dim_factor(number_of_variables) + penalty_factor * penalty + fopt;
  }

private:
  static constexpr double penalty_factor = 1.0;

  const vars_shift_t shift;
  const vars_permblockdiag_t rotation1;
  const vars_conditioning_t conditioning;
  const vars_round_step_t round_step;
  const vars_permblockdiag_t rotation2;
  std::vector<double> scales; /**< @brief condition^(i / (n - 1)) */
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};
