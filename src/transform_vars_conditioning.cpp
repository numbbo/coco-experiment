/**
 * @file transform_vars_conditioning.cpp
 * @brief Implementation of conditioning decision values.
 */

#include <math.h>
#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"

/**
 * @brief Conditions the decision values.
 */
class transform_vars_conditioning_t final : public coco_problem_transformed_t {
public:
  transform_vars_conditioning_t(coco_problem_ptr_t inner, const double alpha)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_conditioning"),
        conditioned_x(inner_problem->get_number_of_variables()), scales(inner_problem->get_number_of_variables()) {
    for (size_t i = 0; i < number_of_variables; ++i)
      scales[i] = pow(alpha, 0.5 * static_cast<double>(i) / (static_cast<double>(number_of_variables) - 1.0));
    if (inner_problem->best_parameter_not_zero()) {
      coco_warning("transform_vars_conditioning(): 'best_parameter' not updated, set to NAN");
      inner_problem->set_best_parameter_to_nan();
    }
  }

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(const double *x, double *y) override {

    for (size_t i = 0; i < number_of_variables; ++i)
      conditioned_x[i] = scales[i] * x[i];
    inner_problem->evaluate_function(conditioned_x.data(), y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the gradient of the transformed function.
   */
  void evaluate_gradient(const double *x, double *y) override {
    for (size_t i = 0; i < number_of_variables; ++i)
      conditioned_x[i] = scales[i] * x[i];
    inner_problem->evaluate_gradient(conditioned_x.data(), y);
    for (size_t i = 0; i < number_of_variables; ++i)
      y[i] = scales[i] * y[i];
  }

  /** @brief Divides by the scaling coefficients. */
  bool inverse(const double *inner_x, double *x) const override {
    for (size_t i = 0; i < number_of_variables; ++i)
      x[i] = inner_x[i] / scales[i];
    return true;
  }

private:
  std::vector<double> conditioned_x; /**< @brief The conditioned decision values */
  std::vector<double> scales;        /**< @brief alpha^(0.5 i / (n - 1)) */
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_conditioning(coco_problem_ptr_t inner_problem, const double alpha) {
  return std::make_unique<transform_vars_conditioning_t>(std::move(inner_problem), alpha);
}
