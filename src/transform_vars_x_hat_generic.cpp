/**
 * @file transform_vars_x_hat_generic.cpp
 * @brief Implementation of multiplying the decision values by the vector 1+-.
 * Wassim: TODO: should eventually replace the non generic version in its use in Schwefel where xopt would be set
 * elsewhere
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

/**
 * @brief Data that the large-scale Lunacek bi-Rastrigin function shares with transform_vars_x_hat_generic.
 */
struct f_lunacek_bi_rastrigin_shared_data_t {
  coco_problem_t *sub_problem_mu0 = nullptr;
  coco_problem_t *sub_problem_mu1 = nullptr;
  std::vector<double> x_hat;
};

/**
 * @brief Multiplies the decision values by 2 and a given vector of 1 and -1.
 *
 * Also stores the result in x_hat of the data shared with the Lunacek bi-Rastrigin function.
 */
struct transform_vars_x_hat_generic_t final : public coco_problem_transformed_t {
  const std::vector<double> sign_vector;
  std::vector<double> x_hat; /**< @brief The transformed decision values */
  f_lunacek_bi_rastrigin_shared_data_t *const shared; /**< @brief Owned by the Lunacek bi-Rastrigin function */

  transform_vars_x_hat_generic_t(coco_problem_t *inner_problem, const double *sign_vector,
                                 f_lunacek_bi_rastrigin_shared_data_t *shared)
      : coco_problem_transformed_t(inner_problem, "transform_vars_x_hat_generic"),
        sign_vector(sign_vector, sign_vector + inner_problem->number_of_variables),
        x_hat(inner_problem->number_of_variables), shared(shared) {}

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(const double *x, double *y) override {
    size_t i;

    for (i = 0; i < number_of_variables; ++i) {
      x_hat[i] = 2.0 * sign_vector[i] * x[i];
      shared->x_hat[i] = x_hat[i];
    }
    coco_evaluate_function(inner_problem, x_hat.data(), y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_vars_x_hat_generic(coco_problem_t *inner_problem, const double *sign_vector,
                                                    f_lunacek_bi_rastrigin_shared_data_t *shared) {
  return new transform_vars_x_hat_generic_t(inner_problem, sign_vector, shared);
}
