/**
 * @file f_bueche_rastrigin.cpp
 * @brief Implementation of the Bueche-Rastrigin function and problem.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

namespace coco {

/**
 * @brief Implements the Bueche-Rastrigin function without connections to any
 * COCO structures.
 */
static double BuecheRastriginRaw(std::span<const double> x) {
  const size_t number_of_variables = x.size();
  double tmp = 0., tmp2 = 0.;

  for (size_t i = 0; i < number_of_variables; ++i) {
    tmp += cos(2 * coco_pi * x[i]);
    tmp2 += x[i] * x[i];
  }
  double result = 10.0 * (static_cast<double>(number_of_variables) - tmp) + tmp2 + 0;
  return result;
}

/**
 * @brief The Bueche-Rastrigin function.
 */
class BuecheRastrigin final : public Problem {
public:
  explicit BuecheRastrigin(const size_t number_of_variables) : Problem(number_of_variables, 1, 0) {
    InitFromScalars("Bueche-Rastrigin function", -5.0, 5.0, 0.0);
    SetId("%s_d%02zu", "bueche-rastrigin", number_of_variables);

    /* Compute best solution */
    EvaluateFunction(best_parameter_, best_value_);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives_ == 1);
    y[0] = BuecheRastriginRaw(x);
    assert(y[0] + 1e-13 >= best_value_[0]);
  }
};

/**
 * @brief The BBOB Bueche-Rastrigin function f4, also in the bbob-boxed and bbob-largescale suites.
 *
 * f(x) = f_rastrigin(s T_osz(x - xopt)) + fopt + 100 f_pen(x) with the 's_i scaling' s, normalized by the
 * dimension before fopt is added in the large-scale suite.
 */
class BbobBuecheRastrigin final : public BbobProblem<BbobBuecheRastrigin> {
public:
  BbobBuecheRastrigin(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                      const BbobVariant variant, const char *problem_id_template, const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "1-separable", problem_id_template, problem_name_template),
        shift_(Xopt(rseed, dimension, variant)), brs_(dimension), fopt_(Bbob2009ComputeFopt(function, instance)),
        normalize_(variant == BbobVariant::LARGESCALE), z_(dimension) {
    /* The optimum of the Bueche-Rastrigin function is the origin, map it back through the transformations */
    brs_.Inverse(best_parameter_, z_);
    oscillate_.Inverse(z_, z_);
    shift_.Inverse(z_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = BoundaryPenalty(x, smallest_values_of_interest_, largest_values_of_interest_);
    shift_.Apply(x, z_);
    oscillate_.Apply(z_, z_);
    brs_.Apply(z_, z_);
    double y = BuecheRastriginRaw(z_);
    if (normalize_)
      y *= NormByDimFactor(number_of_variables_);
    y += fopt_;
    return y + PENALTY_FACTOR * penalty;
  }

private:
  /** @brief The location of the optimum: xopt with the absolute value of every other component. */
  static std::vector<double> Xopt(const long rseed, const size_t dimension, const BbobVariant variant) {
    std::vector<double> xopt = BbobXopt(rseed, dimension, variant);
    /* OME: This step is in the legacy C code but _not_ in the function
     * description. */
    for (size_t i = 0; i < dimension; i += 2) {
      xopt[i] = fabs(xopt[i]);
    }
    return xopt;
  }

  static constexpr double PENALTY_FACTOR = 100.0;

  const VarsShift shift_;
  const VarsOscillate oscillate_{};
  const VarsBrs brs_;
  const double fopt_;
  const bool normalize_;
  std::vector<double> z_; /**< @brief The transformed decision values */
};

} // namespace coco
