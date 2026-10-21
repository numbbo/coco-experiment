/**
 * @file f_schwefel.cpp
 * @brief Implementation of the Schwefel function and problem.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "bbob_problem.cpp"
#include "f_schwefel_generalized.cpp"

/**
 * @brief Implements the Schwefel function without connections to any COCO
 * structures.
 */
static double f_schwefel_raw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  size_t i = 0;
  double penalty, sum;

  /* Boundary handling*/
  penalty = 0.0;
  for (i = 0; i < number_of_variables; ++i) {
    const double tmp = fabs(x[i]) - 500.0;
    if (tmp > 0.0)
      penalty += tmp * tmp;
  }

  /* Computation core */
  sum = 0.0;
  for (i = 0; i < number_of_variables; ++i) {
    sum += x[i] * sin(sqrt(fabs(x[i])));
  }
  double result = 0.01 * (penalty + 418.9828872724339 - sum / static_cast<double>(number_of_variables));

  return result;
}

/**
 * @brief The BBOB Schwefel function f20, also in the bbob-boxed and bbob-largescale suites; the large-scale suite
 * uses f_schwefel_generalized_raw().
 *
 * f(x) = f_schwefel(100 (Λ^10 (ẑ - 2 |xopt|) + 2 |xopt|)) + fopt with ẑ = z^hat(2 x^hat), where x^hat flips the
 * signs of x at random, and xopt = ±4.2096874637 / 2 with the same signs.
 */
class bbob_schwefel_t final : public bbob_problem_t<bbob_schwefel_t> {
public:
  bbob_schwefel_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                  const coco_bbob_variant_t variant, const char *problem_id_template, const char *problem_name_template)
      : bbob_schwefel_t(function, dimension, instance, rseed, variant, schwefel_xopt(rseed, dimension),
                        problem_id_template, problem_name_template) {}

  double operator()(std::span<const double> x) {
    x_hat.apply(x, z1);
    doubling.apply(z1, z1);
    z_hat.apply(z1, z2);
    shift.apply(z2, z2);
    conditioning.apply(z2, z2);
    shift_back.apply(z2, z2);
    scale.apply(z2, z2);
    return raw(z2) + fopt;
  }

private:
  bbob_schwefel_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                  const coco_bbob_variant_t variant, const std::vector<double> &xopt, const char *problem_id_template,
                  const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "5-weakly-structured", problem_id_template,
                       problem_name_template),
        x_hat(rseed, dimension), doubling(2, dimension), z_hat(xopt), shift(two_abs(xopt, 1.0)),
        conditioning(10.0, dimension), shift_back(two_abs(xopt, -1.0)), scale(100, dimension),
        fopt(bbob2009_compute_fopt(function, instance)), generalized(variant == coco_bbob_variant_t::largescale),
        z1(dimension), z2(dimension) {
    /* The best parameter of the Schwefel function (420.96874633 / 100) and the shift (2 |xopt[i]| = 4.2096874637)
     * differ by 4e-10. The best parameter is mapped back through the transformations except the conditioning and
     * the z^hat transformation, which would change it by about 1e-10: the bi-objective suites evaluate their other
     * objective at this point to compute the nadir, which has to stay the same for the existing reference data. At
     * the optimum, the z^hat transformation is the identity. */
    const std::vector<double> optimum(dimension, 420.96874633);
    scale.inverse(optimum, z1);
    shift_back.inverse(z1, z1);
    shift.inverse(z1, z1);
    doubling.inverse(z1, z1);
    x_hat.inverse(z1, best_parameter);
    best_value[0] = raw(optimum) + fopt;
  }

  /** @brief Returns xopt = ±4.2096874637 / 2 with random signs. */
  static std::vector<double> schwefel_xopt(const long rseed, const size_t dimension) {
    std::vector<double> xopt(dimension);
    bbob2009_unif(xopt, rseed);
    for (size_t i = 0; i < dimension; ++i)
      xopt[i] = (xopt[i] < 0.5 ? -1 : 1) * 0.5 * 4.2096874637;
    return xopt;
  }

  /** @brief Returns sign 2 |xopt|. */
  static std::vector<double> two_abs(const std::vector<double> &xopt, const double sign) {
    std::vector<double> result(xopt.size());
    for (size_t i = 0; i < xopt.size(); ++i)
      result[i] = sign * 2 * fabs(xopt[i]);
    return result;
  }

  double raw(std::span<const double> x) const {
    return generalized ? f_schwefel_generalized_raw(x) : f_schwefel_raw(x);
  }

  const vars_x_hat_t x_hat;
  const vars_scale_t doubling;
  const vars_z_hat_t z_hat;
  const vars_shift_t shift; /**< @brief Subtracts 2 |xopt| */
  const vars_conditioning_t conditioning;
  const vars_shift_t shift_back; /**< @brief Adds 2 |xopt| */
  const vars_scale_t scale;
  const double fopt;
  const bool generalized;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};
