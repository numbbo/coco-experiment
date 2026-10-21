/**
 * @file transform_obj_power.cpp
 * @brief Implementation of raising the objective value to the power of a given exponent.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "coco_problem.cpp"

/**
 * @brief Raises the objective value(s) of the inner problem to the power of an exponent.
 */
struct transform_obj_power_t final : public coco_problem_transformed_t {
  const double exponent;

  transform_obj_power_t(coco_problem_t *inner_problem, const double exponent)
      : coco_problem_transformed_t(inner_problem, "transform_obj_power"), exponent(exponent) {
    /* Compute best value */
    evaluate_function(best_parameter.data(), best_value.data());
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

    coco_evaluate_function(inner_problem, x, y);

    for (i = 0; i < number_of_objectives; i++) {
      y[i] = pow(y[i], exponent);
    }
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_obj_power(coco_problem_t *inner_problem, const double exponent) {
  return new transform_obj_power_t(inner_problem, exponent);
}
