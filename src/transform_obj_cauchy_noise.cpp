/**
 * @file transform_obj_cauchy_noise.cpp
 * @brief Implementation of the Cauchy noise model
 */

#include "coco.h"
#include "suite_bbob_noisy_utilities.cpp"

namespace coco {

/**
 * @brief Adds Cauchy noise to the objective value of the inner problem.
 */
class TransformObjCauchyNoise final : public TransformedProblem {
public:
  TransformObjCauchyNoise(ProblemPtr inner, const double alpha, const double p)
      : TransformedProblem(std::move(inner), "cauchy_noise_model"), alpha_(alpha), p_(p) {
    noisy_ = true;
  }

  /**
   * @brief Evaluates the transformed objective function by applying cauchy additive noise.
   */
  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    double uniform_indicator, numerator_normal_variate, denominator_normal_variate;
    uniform_indicator = SampleUniformNoise();
    numerator_normal_variate = SampleGaussianNoise();
    denominator_normal_variate = SampleGaussianNoise();
    denominator_normal_variate = fabs(denominator_normal_variate + 1e-199);
    double cauchy_noise = numerator_normal_variate / (denominator_normal_variate);
    cauchy_noise = uniform_indicator < p_ ? 1e3 + cauchy_noise : 1e3;
    cauchy_noise = alpha_ * cauchy_noise;
    cauchy_noise = cauchy_noise > 0 ? cauchy_noise : 0.;
    double tol = 1e-8;
    inner_problem_->EvaluateFunction(x, y);
    for (size_t i = 0; i < number_of_objectives_; i++) {
      last_noise_free_values_[i] = y[i];
    }
    y[0] = y[0] + cauchy_noise + 1.01 * tol + NoisyBoundaryPenalty(x);
  }

  /** @brief The decision values are not transformed. */
  bool Inverse(std::span<const double> inner_x, std::span<double> x) const override {
    return IdentityInverse(inner_x, x);
  }

private:
  const double alpha_;
  const double p_;
};

} // namespace coco
