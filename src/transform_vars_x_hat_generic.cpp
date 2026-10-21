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
  coco_problem_ptr_t sub_problem_mu0;
  coco_problem_ptr_t sub_problem_mu1;
  std::vector<double> x_hat;
};

/**
 * @brief Multiplies the decision values by 2 and a given vector of 1 and -1.
 *
 * Also stores the result in x_hat of the data shared with the Lunacek bi-Rastrigin function.
 */
class transform_vars_x_hat_generic_t final : public coco_problem_transformed_t {
public:
  transform_vars_x_hat_generic_t(coco_problem_ptr_t inner, std::span<const double> sign_vector,
                                 f_lunacek_bi_rastrigin_shared_data_t *shared)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_x_hat_generic"),
        sign_vector(sign_vector.begin(), sign_vector.end()), x_hat(inner_problem->get_number_of_variables()),
        shared(shared) {
    update_best_parameter();
  }

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {

    for (size_t i = 0; i < number_of_variables; ++i) {
      x_hat[i] = 2.0 * sign_vector[i] * x[i];
      shared->x_hat[i] = x_hat[i];
    }
    inner_problem->evaluate_function(x_hat, y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /** @brief Divides by twice the sign. */
  bool inverse(std::span<const double> inner_x, std::span<double> x) const override {
    for (size_t i = 0; i < number_of_variables; ++i)
      x[i] = inner_x[i] / (2.0 * sign_vector[i]);
    return true;
  }

private:
  const std::vector<double> sign_vector;
  std::vector<double> x_hat;                          /**< @brief The transformed decision values */
  f_lunacek_bi_rastrigin_shared_data_t *const shared; /**< @brief Owned by the Lunacek bi-Rastrigin function */
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_x_hat_generic(coco_problem_ptr_t inner_problem,
                                                       std::span<const double> sign_vector,
                                                       f_lunacek_bi_rastrigin_shared_data_t *shared) {
  return std::make_unique<transform_vars_x_hat_generic_t>(std::move(inner_problem), sign_vector, shared);
}
