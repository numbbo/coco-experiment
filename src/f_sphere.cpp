/**
 * @file f_sphere.cpp
 * @brief Implementation of the sphere function and problem.
 */

#include <assert.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

namespace coco {

/**
 * @brief Implements the sphere function without connections to any COCO
 * structures.
 */
static double SphereRaw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  size_t i = 0;

  double result = 0.0;
  for (i = 0; i < number_of_variables; ++i) {
    result += x[i] * x[i];
  }

  return result;
}

/**
 * @brief The sphere function.
 */
class Sphere final : public Problem {
public:
  explicit Sphere(const size_t number_of_variables) : Problem(number_of_variables, 1, 0) {
    InitFromScalars("sphere function", -5.0, 5.0, 0.0);
    SetId("%s_d%02zu", "sphere", number_of_variables);

    /* Compute best solution */
    EvaluateFunction(best_parameter_, best_value_);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives_ == 1);
    y[0] = SphereRaw(x);
    assert(y[0] + 1e-13 >= best_value_[0]);
  }

  /**
   * @brief Evaluates the gradient of the sphere function.
   */
  void EvaluateGradient(std::span<const double> x, std::span<double> y) override {

    for (size_t i = 0; i < number_of_variables_; ++i) {
      y[i] = 2.0 * x[i];
    }
  }
};

/**
 * @brief The BBOB sphere function f1, also in the bbob-boxed, bbob-largescale and bbob-noisy suites.
 *
 * f(x) = f_sphere(x - xopt) + fopt, normalized by the dimension before fopt is added in the large-scale suite.
 */
class BbobSphere final : public BbobProblem<BbobSphere> {
public:
  BbobSphere(const size_t function, const size_t dimension, const size_t instance, const long rseed,
             const BbobVariant variant, const char *problem_id_template, const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "1-separable", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)), fopt_(Bbob2009ComputeFopt(function, instance)),
        normalize_(variant == BbobVariant::LARGESCALE), z_(dimension) {
    /* The optimum of the sphere function is the origin, map it back through the transformation */
    shift_.Inverse(best_parameter_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    shift_.Apply(x, z_);
    double y = SphereRaw(z_);
    if (normalize_)
      y *= NormByDimFactor(number_of_variables_);
    return y + fopt_;
  }

private:
  const VarsShift shift_;
  const double fopt_;
  const bool normalize_;
  std::vector<double> z_; /**< @brief The transformed decision values */
};

} // namespace coco
