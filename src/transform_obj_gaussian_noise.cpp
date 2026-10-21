/**
 * @file transform_obj_gaussian_noise.cpp
 * @brief Implementation of the Gaussian noise model
 */

#include "coco.h"
#include "suite_bbob_noisy_utilities.cpp"

namespace coco {

/**
 * @brief Applies multiplicative Gaussian noise to the objective value of the inner problem.
 */
class TransformObjGaussianNoise final : public TransformedProblem {
public:
  TransformObjGaussianNoise(ProblemPtr inner, const double beta)
      : TransformedProblem(std::move(inner), "gaussian_noise_model"), beta_(beta) {
    noisy_ = true;
  }

  /**
   * @brief Evaluates the transformed objective function by applying gaussian multiplicative noise.
   */
  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    double fopt = inner_problem_->best_value()[0];
    double gaussian_noise = SampleGaussianNoise();
    gaussian_noise = exp(beta_ * gaussian_noise);
    double tol = 1e-8;
    inner_problem_->EvaluateFunction(x, y);
    for (size_t i = 0; i < number_of_objectives_; i++) {
      last_noise_free_values_[i] = y[i];
    }
    y[0] = y[0] - fopt;
    y[0] = y[0] * gaussian_noise + 1.01 * tol;
    y[0] = y[0] + fopt + NoisyBoundaryPenalty(x);
  }

  /** @brief The decision values are not transformed. */
  bool Inverse(std::span<const double> inner_x, std::span<double> x) const override {
    return IdentityInverse(inner_x, x);
  }

private:
  const double beta_;
};

} // namespace coco
