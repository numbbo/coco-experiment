/**
 * @file f_sharp_ridge_generalized.cpp
 * @brief Implementation of the generalized sharp ridge function and problem.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

namespace coco {

/**
 * @brief The sharp ridge function.
 */
class SharpRidgeGeneralized final : public Problem {
public:
  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double Raw(std::span<const double> x) const;

  SharpRidgeGeneralized(const size_t number_of_variables, size_t proportion_of_linear_dims)
      : Problem(number_of_variables, 1, 0) {
    /* Wassim: proportion_of_linear_dims should probably be allowed to be non-integer */
    InitFromScalars("sharp ridge function", -5.0, 5.0, 0.0);
    SetId("%s_d%02zu", "sharp_ridge_generalized", number_of_variables);
    this->proportion_of_linear_dims_ = proportion_of_linear_dims;

    /* Compute best solution */
    EvaluateFunction(best_parameter_, best_value_);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives_ == 1);
    y[0] = Raw(x);
    assert(y[0] + 1e-13 >= best_value_[0]);
  }

private:
  size_t proportion_of_linear_dims_{};
};

/**
 * @brief Implements the generalized sharp ridge function without connections to
 * any COCO structures.
 */
double SharpRidgeGeneralized::Raw(std::span<const double> x) const {

  static const double alpha = 100.0;
  size_t i = 0, number_linear_dimensions;

  assert(number_of_variables_ > 1);

  double result = 0.0;
  number_linear_dimensions = number_of_variables_ / proportion_of_linear_dims_;
  if (number_of_variables_ % proportion_of_linear_dims_ != 0) {
    number_linear_dimensions += 1;
  }
  for (i = number_linear_dimensions; i < number_of_variables_; ++i) {
    result += x[i] * x[i];
  }
  result = alpha * sqrt(result);
  for (i = 0; i < number_linear_dimensions; i++) {
    result += x[i] * x[i];
  }

  return result;
}

/**
 * @brief The sharp ridge function f13 of the bbob-largescale suite, generalized to n / 40 linear dimensions and
 * with permuted block rotations.
 *
 * f(x) = f_sharp_ridge_generalized(Q Λ^10 R (x - xopt)) / n + fopt with the permuted block rotations R and Q.
 */
class BbobSharpRidgeGeneralizedPermblockdiag final : public BbobProblem<BbobSharpRidgeGeneralizedPermblockdiag> {
public:
  BbobSharpRidgeGeneralizedPermblockdiag(const size_t function, const size_t dimension, const size_t instance,
                                         const long rseed, const BbobVariant variant, const char *problem_id_template,
                                         const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "3-ill-conditioned", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)), rotation1_(dimension, rseed + 5000000, rseed, rseed + 6000000),
        conditioning_(10.0, dimension), rotation2_(dimension, rseed + 3000000, rseed + 1000000, rseed + 4000000),
        sharp_ridge_(dimension, 40), fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension), z2_(dimension) {
    /* The optimum of the sharp ridge function is the origin, map it back through the transformations */
    rotation2_.Inverse(best_parameter_, z1_);
    conditioning_.Inverse(z1_, z1_);
    rotation1_.Inverse(z1_, z2_);
    shift_.Inverse(z2_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    shift_.Apply(x, z1_);
    rotation1_.Apply(z1_, z2_);
    conditioning_.Apply(z2_, z2_);
    rotation2_.Apply(z2_, z1_);
    return sharp_ridge_.Raw(z1_) * NormByDimFactor(number_of_variables_) + fopt_;
  }

private:
  const VarsShift shift_;
  const VarsPermblockdiag rotation1_;
  const VarsConditioning conditioning_;
  const VarsPermblockdiag rotation2_;
  const SharpRidgeGeneralized sharp_ridge_;
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

} // namespace coco
