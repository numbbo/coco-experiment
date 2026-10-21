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

namespace coco {

/**
 * @brief Implements the Schwefel function without connections to any COCO
 * structures.
 */
static double SchwefelRaw(std::span<const double> x) {
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
 * uses SchwefelGeneralizedRaw().
 *
 * f(x) = f_schwefel(100 (Λ^10 (ẑ - 2 |xopt|) + 2 |xopt|)) + fopt with ẑ = z^hat(2 x^hat), where x^hat flips the
 * signs of x at random, and xopt = ±4.2096874637 / 2 with the same signs.
 */
class BbobSchwefel final : public BbobProblem<BbobSchwefel> {
public:
  BbobSchwefel(const size_t function, const size_t dimension, const size_t instance, const long rseed,
               const BbobVariant variant, const char *problem_id_template, const char *problem_name_template)
      : BbobSchwefel(function, dimension, instance, rseed, variant, SchwefelXopt(rseed, dimension), problem_id_template,
                     problem_name_template) {}

  double operator()(std::span<const double> x) {
    x_hat_.Apply(x, z1_);
    doubling_.Apply(z1_, z1_);
    z_hat_.Apply(z1_, z2_);
    shift_.Apply(z2_, z2_);
    conditioning_.Apply(z2_, z2_);
    shift_back_.Apply(z2_, z2_);
    scale_.Apply(z2_, z2_);
    return Raw(z2_) + fopt_;
  }

private:
  BbobSchwefel(const size_t function, const size_t dimension, const size_t instance, const long rseed,
               const BbobVariant variant, const std::vector<double> &xopt, const char *problem_id_template,
               const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "5-weakly-structured", problem_id_template, problem_name_template),
        x_hat_(rseed, dimension), doubling_(2, dimension), z_hat_(xopt), shift_(TwoAbs(xopt, 1.0)),
        conditioning_(10.0, dimension), shift_back_(TwoAbs(xopt, -1.0)), scale_(100, dimension),
        fopt_(Bbob2009ComputeFopt(function, instance)), generalized_(variant == BbobVariant::LARGESCALE),
        z1_(dimension), z2_(dimension) {
    /* The best parameter of the Schwefel function (420.96874633 / 100) and the shift (2 |xopt[i]| = 4.2096874637)
     * differ by 4e-10. The best parameter is mapped back through the transformations except the conditioning and
     * the z^hat transformation, which would change it by about 1e-10: the bi-objective suites evaluate their other
     * objective at this point to compute the nadir, which has to stay the same for the existing reference data. At
     * the optimum, the z^hat transformation is the identity. */
    const std::vector<double> optimum(dimension, 420.96874633);
    scale_.Inverse(optimum, z1_);
    shift_back_.Inverse(z1_, z1_);
    shift_.Inverse(z1_, z1_);
    doubling_.Inverse(z1_, z1_);
    x_hat_.Inverse(z1_, best_parameter_);
    best_value_[0] = Raw(optimum) + fopt_;
  }

  /** @brief Returns xopt = ±4.2096874637 / 2 with random signs. */
  static std::vector<double> SchwefelXopt(const long rseed, const size_t dimension) {
    std::vector<double> xopt(dimension);
    Bbob2009Unif(xopt, rseed);
    for (size_t i = 0; i < dimension; ++i)
      xopt[i] = (xopt[i] < 0.5 ? -1 : 1) * 0.5 * 4.2096874637;
    return xopt;
  }

  /** @brief Returns sign 2 |xopt|. */
  static std::vector<double> TwoAbs(const std::vector<double> &xopt, const double sign) {
    std::vector<double> result(xopt.size());
    for (size_t i = 0; i < xopt.size(); ++i)
      result[i] = sign * 2 * fabs(xopt[i]);
    return result;
  }

  double Raw(std::span<const double> x) const { return generalized_ ? SchwefelGeneralizedRaw(x) : SchwefelRaw(x); }

  const VarsXHat x_hat_;
  const VarsScale doubling_;
  const VarsZHat z_hat_;
  const VarsShift shift_; /**< @brief Subtracts 2 |xopt| */
  const VarsConditioning conditioning_;
  const VarsShift shift_back_; /**< @brief Adds 2 |xopt| */
  const VarsScale scale_;
  const double fopt_;
  const bool generalized_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

} // namespace coco
