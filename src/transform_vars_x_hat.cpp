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
struct transform_vars_x_hat_t final : public coco_problem_transformed_t {
  const long seed;
  std::vector<double> x_hat; /**< @brief The transformed decision values */

  transform_vars_x_hat_t(coco_problem_ptr_t inner, const long seed)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_x_hat"), seed(seed),
        x_hat(inner_problem->number_of_variables) {
    size_t i;
    if (coco_problem_best_parameter_not_zero(this)) {
      bbob2009_unif(x_hat.data(), number_of_variables, seed);
      for (i = 0; i < number_of_variables; ++i)
        if (x_hat[i] < 0.5) /* with probability 1/2 */
          best_parameter[i] *= -1;
    }
  }

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(const double *x, double *y) override {
    size_t i;

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }

    bbob2009_unif(x_hat.data(), number_of_variables, seed);

    for (i = 0; i < number_of_variables; ++i) {
      if (x_hat[i] < 0.5) {
        x_hat[i] = -x[i];
      } else {
        x_hat[i] = x[i];
      }
    }
    coco_evaluate_function(inner_problem.get(), x_hat.data(), y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_x_hat(coco_problem_ptr_t inner_problem, const long seed) {
  return std::make_unique<transform_vars_x_hat_t>(std::move(inner_problem), seed);
}
