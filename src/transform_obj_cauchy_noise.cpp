/**
 * @file transform_obj_cauchy_noise.cpp
 * @brief Implementation of the Cauchy noise model
 */

#include "coco.h"
#include "suite_bbob_noisy_utilities.cpp"

/**
 * @brief Adds Cauchy noise to the objective value of the inner problem.
 */
struct transform_obj_cauchy_noise_t final : public coco_problem_transformed_t {
  const double alpha;
  const double p;

  transform_obj_cauchy_noise_t(coco_problem_t *inner_problem, const double alpha, const double p)
      : coco_problem_transformed_t(inner_problem, "cauchy_noise_model"), alpha(alpha), p(p) {
    is_noisy = 1;
  }

  /**
   * @brief Evaluates the transformed objective function by applying cauchy additive noise.
   */
  void evaluate_function(const double *x, double *y) override {
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
    *(y) = *(y) + cauchy_noise + 1.01 * tol + coco_boundary_handling(this, x);
  }
};

/**
 * @brief Allocates a noisy problem with cauchy noise.
 */
static coco_problem_t *transform_obj_cauchy_noise(coco_problem_t *inner_problem, const double alpha, const double p) {
  return new transform_obj_cauchy_noise_t(inner_problem, alpha, p);
}
