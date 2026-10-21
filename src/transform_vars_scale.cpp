/**
 * @file transform_vars_scale.cpp
 * @brief Implementation of scaling decision values by a given factor.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"

/**
 * @brief Scales the decision values by a factor.
 */
struct transform_vars_scale_t final : public coco_problem_transformed_t {
  const double factor;
  std::vector<double> scaled_x; /**< @brief The scaled decision values */

  transform_vars_scale_t(coco_problem_t *inner_problem, const double factor)
      : coco_problem_transformed_t(inner_problem, "transform_vars_scale"), factor(factor),
        scaled_x(inner_problem->number_of_variables) {
    size_t i;
    /* Compute best parameter */
    if (factor != 0.) {
      for (i = 0; i < number_of_variables; i++) {
        best_parameter[i] /= factor;
      }
    } /* else error? */
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
      scaled_x[i] = factor * x[i];
    }
    coco_evaluate_function(inner_problem, scaled_x.data(), y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_vars_scale(coco_problem_t *inner_problem, const double factor) {
  return new transform_vars_scale_t(inner_problem, factor);
}
