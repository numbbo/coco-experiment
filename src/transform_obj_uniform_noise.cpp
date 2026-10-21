/**
 * @file transform_obj_uniform_noise.cpp
 * @brief Implementation of the Uniform noise model
 */

#include "coco.h"
#include "suite_bbob_noisy_utilities.cpp"

namespace coco {

/**
 * @brief Applies multiplicative uniform noise to the objective value of the inner problem.
 */
class TransformObjUniformNoise final : public TransformedProblem {
public:
  TransformObjUniformNoise(ProblemPtr inner, const double alpha, const double beta)
      : TransformedProblem(std::move(inner), "uniform_noise_model"), alpha_(alpha), beta_(beta) {
    noisy_ = true;
  }

  /**
   * @brief Evaluates the transformed objective function by applying uniform multiplicative noise.
   */
  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    double uniform_noise_term1, uniform_noise_term2;
    double fopt = inner_problem_->best_value()[0];
    uniform_noise_term1 = SampleUniformNoise();
    uniform_noise_term2 = SampleUniformNoise();
    double uniform_noise_factor = pow(uniform_noise_term1, beta_);
    inner_problem_->EvaluateFunction(x, y);
    for (size_t i = 0; i < number_of_objectives_; i++) {
      last_noise_free_values_[i] = y[i];
    }
    y[0] = y[0] - fopt;
    double scaling_factor = 1e9 / (y[0] + 1e-99);
    scaling_factor = pow(scaling_factor, alpha_ * uniform_noise_term2);
    scaling_factor = scaling_factor > 1 ? scaling_factor : 1;
    double uniform_noise = uniform_noise_factor * scaling_factor;
    double tol = 1e-8;
    y[0] = y[0] * uniform_noise + 1.01 * tol;
    y[0] = y[0] + fopt + NoisyBoundaryPenalty(x);
  }

  /** @brief The decision values are not transformed. */
  bool Inverse(std::span<const double> inner_x, std::span<double> x) const override {
    return IdentityInverse(inner_x, x);
  }

private:
  const double alpha_;
  const double beta_;
};

} // namespace coco
