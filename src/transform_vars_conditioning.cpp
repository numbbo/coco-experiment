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
struct transform_vars_conditioning_t final : public coco_problem_transformed_t {
  std::vector<double> conditioned_x; /**< @brief The conditioned decision values */
  const double alpha;

  transform_vars_conditioning_t(coco_problem_ptr_t inner, const double alpha)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_conditioning"),
        conditioned_x(inner_problem->number_of_variables), alpha(alpha) {
    if (inner_problem->best_parameter_not_zero()) {
      coco_warning("transform_vars_conditioning(): 'best_parameter' not updated, set to NAN");
      coco_vector_set_to_nan(inner_problem->best_parameter.data(), inner_problem->number_of_variables);
    }
  }

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(const double *x, double *y) override {

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }

    for (size_t i = 0; i < number_of_variables; ++i) {
      /* OME: We could precalculate the scaling coefficients if we
       * really wanted to.
       */
      conditioned_x[i] = pow(alpha, 0.5 * (double)(long)i / ((double)(long)number_of_variables - 1.0)) * x[i];
    }
    coco_evaluate_function(inner_problem.get(), conditioned_x.data(), y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the gradient of the transformed function.
   */
  void evaluate_gradient(const double *x, double *y) override {

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }

    std::vector<double> gradient(inner_problem->number_of_variables, 0.0);

    for (size_t i = 0; i < number_of_variables; ++i) {
      gradient[i] = pow(alpha, 0.5 * (double)(long)i / ((double)(long)number_of_variables - 1.0));
      conditioned_x[i] = gradient[i] * x[i];
    }
    inner_problem->evaluate_gradient(conditioned_x.data(), y);

    for (size_t i = 0; i < inner_problem->number_of_variables; ++i)
      gradient[i] *= y[i];

    for (size_t i = 0; i < inner_problem->number_of_variables; ++i)
      y[i] = gradient[i];

  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_conditioning(coco_problem_ptr_t inner_problem, const double alpha) {
  return std::make_unique<transform_vars_conditioning_t>(std::move(inner_problem), alpha);
}
