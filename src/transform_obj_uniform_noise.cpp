/**
 * @file transform_obj_uniform_noise.cpp
 * @brief Implementation of the Uniform noise model
 */

#include "coco.h"
#include "suite_bbob_noisy_utilities.cpp"

/**
 * @brief Applies multiplicative uniform noise to the objective value of the inner problem.
 */
class transform_obj_uniform_noise_t final : public coco_problem_transformed_t {
public:
  transform_obj_uniform_noise_t(coco_problem_ptr_t inner, const double alpha, const double beta)
      : coco_problem_transformed_t(std::move(inner), "uniform_noise_model"), alpha(alpha), beta(beta) {
    noisy = true;
  }

  /**
   * @brief Evaluates the transformed objective function by applying uniform multiplicative noise.
   */
  void evaluate_function(const double *x, double *y) override {
    double uniform_noise_term1, uniform_noise_term2;
    double fopt = inner_problem->get_best_value()[0];
    uniform_noise_term1 = coco_sample_uniform_noise();
    uniform_noise_term2 = coco_sample_uniform_noise();
    double uniform_noise_factor = pow(uniform_noise_term1, beta);
    inner_problem->evaluate_function(x, y);
    for (size_t i = 0; i < number_of_objectives; i++) {
      last_noise_free_values[i] = y[i];
    }
    *(y) = *(y)-fopt;
    double scaling_factor = 1e9 / (*(y) + 1e-99);
    scaling_factor = pow(scaling_factor, alpha * uniform_noise_term2);
    scaling_factor = scaling_factor > 1 ? scaling_factor : 1;
    double uniform_noise = uniform_noise_factor * scaling_factor;
    double tol = 1e-8;
    *(y) = *(y)*uniform_noise + 1.01 * tol;
    *(y) = *(y) + fopt + boundary_penalty(x);
  }

  /** @brief The decision values are not transformed. */
  bool inverse(const double *inner_x, double *x) const override { return identity_inverse(inner_x, x); }

private:
  const double alpha;
  const double beta;
};

/**
 * @brief Allocates a noisy problem with uniform noise.
 */
static coco_problem_ptr_t transform_obj_uniform_noise(coco_problem_ptr_t inner_problem, const double alpha,
                                                      const double beta) {
  return std::make_unique<transform_obj_uniform_noise_t>(std::move(inner_problem), alpha, beta);
}
