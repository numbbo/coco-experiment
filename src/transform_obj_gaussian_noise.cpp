/**
 * @file transform_obj_gaussian_noise.cpp
 * @brief Implementation of the Gaussian noise model
 */

#include "coco.h"
#include "suite_bbob_noisy_utilities.cpp"

/**
 * @brief Applies multiplicative Gaussian noise to the objective value of the inner problem.
 */
class transform_obj_gaussian_noise_t final : public coco_problem_transformed_t {
public:
  transform_obj_gaussian_noise_t(coco_problem_ptr_t inner, const double beta)
      : coco_problem_transformed_t(std::move(inner), "gaussian_noise_model"), beta(beta) {
    noisy = true;
  }

  /**
   * @brief Evaluates the transformed objective function by applying gaussian multiplicative noise.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    double fopt = inner_problem->get_best_value()[0];
    double gaussian_noise = coco_sample_gaussian_noise();
    gaussian_noise = exp(beta * gaussian_noise);
    double tol = 1e-8;
    inner_problem->evaluate_function(x, y);
    for (size_t i = 0; i < number_of_objectives; i++) {
      last_noise_free_values[i] = y[i];
    }
    y[0] = y[0] - fopt;
    y[0] = y[0] * gaussian_noise + 1.01 * tol;
    y[0] = y[0] + fopt + boundary_penalty(x);
  }

  /** @brief The decision values are not transformed. */
  bool inverse(std::span<const double> inner_x, std::span<double> x) const override {
    return identity_inverse(inner_x, x);
  }

private:
  const double beta;
};

/**
 * @brief Allocates a noisy problem with gaussian noise.
 */
static coco_problem_ptr_t transform_obj_gaussian_noise(coco_problem_ptr_t inner_problem, const double beta) {
  return std::make_unique<transform_obj_gaussian_noise_t>(std::move(inner_problem), beta);
}
