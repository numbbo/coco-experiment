/**
 * @file f_lunacek_bi_rastrigin.cpp
 * @brief Implementation of the Lunacek bi-Rastrigin function and problem.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "bbob_problem.cpp"
#include "f_sphere.cpp"

/**
 * @brief The BBOB Lunacek bi-Rastrigin function f24, also in the bbob-boxed suite.
 *
 * f(x) = min(Σ (x^_i - μ0)^2, d n + s Σ (x^_i - μ1)^2) + 10 (n - Σ cos(2π z_i)) + 10^4 f_pen(x) + fopt with
 * x^ = 2 sign(xopt) x and z = Q Λ^100 R (x^ - μ0).
 */
class bbob_lunacek_bi_rastrigin_t final : public bbob_problem_t<bbob_lunacek_bi_rastrigin_t> {
public:
  bbob_lunacek_bi_rastrigin_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                              const coco_bbob_variant_t variant, const char *problem_id_template,
                              const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "5-weakly-structured", problem_id_template,
                       problem_name_template),
        xopt(lunacek_xopt(rseed, dimension)), shift(std::vector<double>(dimension, mu0)),
        conditioned_rotation(conditioned_rotation_matrix(rseed, dimension), std::vector<double>(dimension)),
        rotation(bbob2009_compute_rotation(rseed + 1000000, dimension), std::vector<double>(dimension)),
        s(1. - 0.5 / (sqrt(static_cast<double>(dimension + 20)) - 4.1)), mu1(-sqrt((mu0 * mu0 - d) / s)),
        fopt(bbob2009_compute_fopt(function, instance)), x_hat(dimension), z1(dimension), z2(dimension) {
    /* The variant only matters for the xopt that the legacy code computes and discards */
    (void)variant;
    best_parameter = xopt;
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = coco_boundary_penalty(x, smallest_values_of_interest, largest_values_of_interest);
    double sum1 = 0., sum2 = 0., sum3 = 0.;

    for (size_t i = 0; i < number_of_variables; ++i) {
      x_hat[i] = 2. * x[i];
      if (xopt[i] < 0.) {
        x_hat[i] *= -1.;
      }
    }
    shift.apply(x_hat, z1);
    conditioned_rotation.apply(z1, z2);
    rotation.apply(z2, z1);
    for (size_t i = 0; i < number_of_variables; ++i) {
      sum1 += (x_hat[i] - mu0) * (x_hat[i] - mu0);
      sum2 += (x_hat[i] - mu1) * (x_hat[i] - mu1);
      sum3 += cos(2 * coco_pi * z1[i]);
    }
    const double y = coco_double_min(sum1, d * static_cast<double>(number_of_variables) + s * sum2) +
                     10. * (static_cast<double>(number_of_variables) - sum3) + 1e4 * penalty;
    return y + fopt;
  }

private:
  static constexpr double mu0 = 2.5;
  static constexpr double d = 1.;

  /** @brief Returns xopt = ±μ0 / 2 with random signs. */
  static std::vector<double> lunacek_xopt(const long rseed, const size_t dimension) {
    std::vector<double> signs(dimension), xopt(dimension);
    bbob2009_gauss(signs, rseed);
    for (size_t i = 0; i < dimension; ++i) {
      xopt[i] = 0.5 * mu0;
      if (signs[i] < 0.0) {
        xopt[i] *= -1.0;
      }
    }
    return xopt;
  }

  /** @brief Returns Λ^100 R with the rotation R. */
  static coco_matrix_t conditioned_rotation_matrix(const long rseed, const size_t dimension) {
    static const double condition = 100.;
    assert(dimension > 1);
    coco_matrix_t M = bbob2009_compute_rotation(rseed, dimension);
    for (size_t i = 0; i < dimension; ++i) {
      const double c1 = pow(sqrt(condition), (static_cast<double>(i)) / static_cast<double>(dimension - 1));
      for (size_t j = 0; j < dimension; ++j)
        M(i, j) = c1 * M(i, j);
    }
    return M;
  }

  const std::vector<double> xopt;
  const vars_shift_t shift; /**< @brief Subtracts μ0 */
  const vars_affine_t conditioned_rotation;
  const vars_affine_t rotation;
  const double s, mu1;
  const double fopt;
  std::vector<double> x_hat, z1, z2; /**< @brief The transformed decision values */
};

/**
 * @brief The Lunacek bi-Rastrigin function f24 of the bbob-largescale suite, with permuted block rotations.
 *
 * f(x) = (10 (n - Σ cos(2π z_i)) + min(f_sphere(x^ - μ0), d n + s f_sphere(x^ - μ1))) / n + 10^4 f_pen(x) + fopt
 * with x^ = 2 sign x, the random signs sign, z = Q Λ^100 R (x^ - μ0) and the permuted block rotations R and Q.
 */
class bbob_lunacek_bi_rastrigin_permblockdiag_t final
    : public bbob_problem_t<bbob_lunacek_bi_rastrigin_permblockdiag_t> {
public:
  bbob_lunacek_bi_rastrigin_permblockdiag_t(const size_t function, const size_t dimension, const size_t instance,
                                            const long rseed, const char *problem_id_template,
                                            const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "5-weakly-structured", problem_id_template,
                       problem_name_template),
        signs(random_signs(rseed, dimension)), s(1. - 0.5 / (sqrt(static_cast<double>(dimension + 20)) - 4.1)),
        mu1(-sqrt((mu0 * mu0 - d) / s)), shift_mu0(std::vector<double>(dimension, mu0)),
        shift_mu1(std::vector<double>(dimension, mu1)), rotation1(dimension, rseed + 3000000, rseed, rseed + 4000000),
        conditioning(100.0, dimension), rotation2(dimension, rseed + 5000000, rseed + 1000000, rseed + 6000000),
        fopt(bbob2009_compute_fopt(function, instance)), x_hat(dimension), z1(dimension), z2(dimension) {
    /* The optimum is where x^ = μ0, map it back through the transformations */
    rotation2.inverse(best_parameter, z1);
    conditioning.inverse(z1, z1);
    rotation1.inverse(z1, z2);
    shift_mu0.inverse(z2, z1);
    for (size_t i = 0; i < dimension; ++i)
      best_parameter[i] = z1[i] / (2.0 * signs[i]);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = coco_boundary_penalty(x, smallest_values_of_interest, largest_values_of_interest);
    for (size_t i = 0; i < number_of_variables; ++i)
      x_hat[i] = 2.0 * signs[i] * x[i];

    shift_mu0.apply(x_hat, z1);
    const double y0 = f_sphere_raw(z1);
    rotation1.apply(z1, z2);
    conditioning.apply(z2, z2);
    rotation2.apply(z2, z1);
    shift_mu1.apply(x_hat, z2);
    const double y1 = f_sphere_raw(z2) * s + d * static_cast<double>(number_of_variables);

    double y = static_cast<double>(number_of_variables);
    for (size_t i = 0; i < number_of_variables; i++) {
      y -= cos(2 * coco_pi * z1[i]);
    }
    y *= 10.0;
    y += coco_double_min(y0, y1);
    y *= coco_norm_by_dim_factor(number_of_variables);
    y += penalty_factor * penalty;
    return y + fopt;
  }

private:
  static constexpr double mu0 = 2.5;
  static constexpr double d = 1.0;
  static constexpr double penalty_factor = 1e4;

  /** @brief Returns ±1 with random signs. */
  static std::vector<double> random_signs(const long rseed, const size_t dimension) {
    std::vector<double> signs(dimension);
    bbob2009_gauss(signs, rseed);
    for (size_t i = 0; i < dimension; i++) {
      signs[i] = signs[i] < 0.0 ? -1.0 : 1.0;
    }
    return signs;
  }

  const std::vector<double> signs;
  const double s, mu1;
  const vars_shift_t shift_mu0; /**< @brief Subtracts μ0 */
  const vars_shift_t shift_mu1; /**< @brief Subtracts μ1 */
  const vars_permblockdiag_t rotation1;
  const vars_conditioning_t conditioning;
  const vars_permblockdiag_t rotation2;
  const double fopt;
  std::vector<double> x_hat, z1, z2; /**< @brief The transformed decision values */
};
