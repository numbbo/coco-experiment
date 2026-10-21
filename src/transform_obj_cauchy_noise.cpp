/**
 * @file transform_obj_cauchy_noise.cpp
 * @brief Implementation of the Cauchy noise model
 */

#include "coco.h"
#include "suite_bbob_noisy_utilities.cpp"

/**
 * @brief Adds Cauchy noise to the objective value of the inner problem.
 */
class transform_obj_cauchy_noise_t final : public coco_problem_transformed_t {
public:
  transform_obj_cauchy_noise_t(coco_problem_ptr_t inner, const double alpha, const double p)
      : coco_problem_transformed_t(std::move(inner), "cauchy_noise_model"), alpha(alpha), p(p) {
    noisy = true;
  }

  /**
   * @brief Evaluates the transformed objective function by applying cauchy additive noise.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    double uniform_indicator, numerator_normal_variate, denominator_normal_variate;
    uniform_indicator = coco_sample_uniform_noise();
    numerator_normal_variate = coco_sample_gaussian_noise();
    denominator_normal_variate = coco_sample_gaussian_noise();
    denominator_normal_variate = fabs(denominator_normal_variate + 1e-199);
    double cauchy_noise = numerator_normal_variate / (denominator_normal_variate);
    cauchy_noise = uniform_indicator < p ? 1e3 + cauchy_noise : 1e3;
    cauchy_noise = alpha * cauchy_noise;
    cauchy_noise = cauchy_noise > 0 ? cauchy_noise : 0.;
    double tol = 1e-8;
    inner_problem->evaluate_function(x, y);
    for (size_t i = 0; i < number_of_objectives; i++) {
      last_noise_free_values[i] = y[i];
    }
    y[0] = y[0] + cauchy_noise + 1.01 * tol + boundary_penalty(x);
  }

  /** @brief The decision values are not transformed. */
  bool inverse(std::span<const double> inner_x, std::span<double> x) const override {
    return identity_inverse(inner_x, x);
  }

private:
  const double alpha;
  const double p;
};

/**
 * @brief Allocates a noisy problem with cauchy noise.
 */
static coco_problem_ptr_t transform_obj_cauchy_noise(coco_problem_ptr_t inner_problem, const double alpha,
                                                     const double p) {
  return std::make_unique<transform_obj_cauchy_noise_t>(std::move(inner_problem), alpha, p);
}
