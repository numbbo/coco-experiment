/**
 * @file f_discus.cpp
 * @brief Implementation of the discus function and problem.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

namespace coco {

/**
 * @brief Implements the discus function without connections to any COCO
 * structures.
 */
static double DiscusRaw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  static const double condition = 1.0e6;

  double result = condition * x[0] * x[0];
  for (size_t i = 1; i < number_of_variables; ++i) {
    result += x[i] * x[i];
  }

  return result;
}

/**
 * @brief The discus function.
 */
class Discus final : public Problem {
public:
  explicit Discus(const size_t number_of_variables) : Problem(number_of_variables, 1, 0) {
    InitFromScalars("discus function", -5.0, 5.0, 0.0);
    SetId("%s_d%02zu", "discus", number_of_variables);

    /* Compute best solution */
    EvaluateFunction(best_parameter_, best_value_);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives_ == 1);
    y[0] = DiscusRaw(x);
    assert(y[0] + 1e-13 >= best_value_[0]);
  }

  /**
   * @brief Evaluates the gradient of the discus function.
   */
  void EvaluateGradient(std::span<const double> x, std::span<double> y) override {

    static const double condition = 1.0e6;

    y[0] = condition * 2.0 * x[0];
    for (size_t i = 1; i < number_of_variables_; ++i)
      y[i] = 2.0 * x[i];
  }
};

/**
 * @brief The BBOB discus function f11, also in the bbob-boxed suite.
 *
 * f(x) = f_discus(T_osz(R (x - xopt))) + fopt.
 */
class BbobDiscus final : public BbobProblem<BbobDiscus> {
public:
  BbobDiscus(const size_t function, const size_t dimension, const size_t instance, const long rseed,
             const BbobVariant variant, const char *problem_id_template, const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "3-ill-conditioned", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)),
        rotation_(Bbob2009ComputeRotation(rseed + 1000000, dimension), std::vector<double>(dimension)),
        fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension), z2_(dimension) {
    /* The optimum of the discus function is the origin, map it back through the transformations */
    oscillate_.Inverse(best_parameter_, z1_);
    rotation_.Inverse(z1_, z2_);
    shift_.Inverse(z2_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    shift_.Apply(x, z1_);
    rotation_.Apply(z1_, z2_);
    oscillate_.Apply(z2_, z2_);
    return DiscusRaw(z2_) + fopt_;
  }

private:
  const VarsShift shift_;
  const VarsAffine rotation_;
  const VarsOscillate oscillate_{};
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

} // namespace coco
