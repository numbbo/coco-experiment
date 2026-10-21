/**
 * @file transform_vars_x_hat.cpp
 * @brief Implementation of multiplying the decision values by the vector 1+-.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

/**
 * @brief Multiplies the decision values by a random vector of 1 and -1.
 */
class transform_vars_x_hat_t final : public coco_problem_transformed_t {
public:
  transform_vars_x_hat_t(coco_problem_ptr_t inner, const long seed)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_x_hat"), seed(seed),
        x_hat(inner_problem->get_number_of_variables()) {
    update_best_parameter();
  }

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {

    bbob2009_unif(x_hat, seed);

    for (size_t i = 0; i < number_of_variables; ++i) {
      if (x_hat[i] < 0.5) {
        x_hat[i] = -x[i];
      } else {
        x_hat[i] = x[i];
      }
    }
    inner_problem->evaluate_function(x_hat, y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /** @brief Flips the same signs again. */
  bool inverse(std::span<const double> inner_x, std::span<double> x) const override {
    std::vector<double> random_numbers(number_of_variables);
    bbob2009_unif(random_numbers, seed);
    for (size_t i = 0; i < number_of_variables; ++i)
      x[i] = random_numbers[i] < 0.5 ? -inner_x[i] : inner_x[i];
    return true;
  }

private:
  const long seed;
  std::vector<double> x_hat; /**< @brief The transformed decision values */
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_x_hat(coco_problem_ptr_t inner_problem, const long seed) {
  return std::make_unique<transform_vars_x_hat_t>(std::move(inner_problem), seed);
}
