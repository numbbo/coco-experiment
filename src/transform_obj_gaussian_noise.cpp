/**
 * @file transform_obj_gaussian_noise.cpp
 * @brief Implementation of the Gaussian noise model
 */

#include "coco.h"
#include "suite_bbob_noisy_utilities.cpp"

/**
 * @brief Applies multiplicative Gaussian noise to the objective value of the inner problem.
 */
struct transform_obj_gaussian_noise_t final : public coco_problem_transformed_t {
  const double beta;

  transform_obj_gaussian_noise_t(coco_problem_t *inner_problem, const double beta)
      : coco_problem_transformed_t(inner_problem, "gaussian_noise_model"), beta(beta) {
    is_noisy = 1;
  }

  /**
   * @brief Evaluates the transformed objective function by applying gaussian multiplicative noise.
   */
  void evaluate_function(const double *x, double *y) override {
    double fopt = inner_problem->best_value[0];
    double gaussian_noise = coco_sample_gaussian_noise();
    gaussian_noise = exp(beta * gaussian_noise);
    double tol = 1e-8;
    inner_problem->evaluate_function(x, y);
    for (size_t i = 0; i < number_of_objectives; i++) {
      last_noise_free_values[i] = y[i];
    }
    *(y) = *(y)-fopt;
    *(y) = *(y)*gaussian_noise + 1.01 * tol;
    *(y) = *(y) + fopt + coco_boundary_handling(this, x);
  }
};

/**
 * @brief Allocates a noisy problem with gaussian noise.
 */
static coco_problem_t *transform_obj_gaussian_noise(coco_problem_t *inner_problem, const double beta) {
  return new transform_obj_gaussian_noise_t(inner_problem, beta);
}
