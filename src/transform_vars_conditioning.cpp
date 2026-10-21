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

  transform_vars_conditioning_t(coco_problem_t *inner_problem, const double alpha)
      : coco_problem_transformed_t(inner_problem, "transform_vars_conditioning"),
        conditioned_x(inner_problem->number_of_variables), alpha(alpha) {
    if (coco_problem_best_parameter_not_zero(inner_problem)) {
      coco_warning("transform_vars_conditioning(): 'best_parameter' not updated, set to NAN");
      coco_vector_set_to_nan(inner_problem->best_parameter.data(), inner_problem->number_of_variables);
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

    for (i = 0; i < number_of_variables; ++i) {
      /* OME: We could precalculate the scaling coefficients if we
       * really wanted to.
       */
      conditioned_x[i] = pow(alpha, 0.5 * (double)(long)i / ((double)(long)number_of_variables - 1.0)) * x[i];
    }
    coco_evaluate_function(inner_problem, conditioned_x.data(), y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the gradient of the transformed function.
   */
  void evaluate_gradient(const double *x, double *y) override {
    size_t i;
    double *gradient;

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }

    gradient = coco_allocate_vector(inner_problem->number_of_variables);

    for (i = 0; i < number_of_variables; ++i) {
      gradient[i] = pow(alpha, 0.5 * (double)(long)i / ((double)(long)number_of_variables - 1.0));
      conditioned_x[i] = gradient[i] * x[i];
    }
    bbob_evaluate_gradient(inner_problem, conditioned_x.data(), y);

    for (i = 0; i < inner_problem->number_of_variables; ++i)
      gradient[i] *= y[i];

    for (i = 0; i < inner_problem->number_of_variables; ++i)
      y[i] = gradient[i];

    coco_free_memory(gradient);
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_vars_conditioning(coco_problem_t *inner_problem, const double alpha) {
  return new transform_vars_conditioning_t(inner_problem, alpha);
}
