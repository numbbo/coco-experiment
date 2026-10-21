/**
 * @file transform_obj_penalize.cpp
 * @brief Implementation of adding a penalty to the objective value for solutions outside of the ROI in the
 * decision space.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"

/**
 * @brief Adds a penalty for solutions outside of the region of interest to the objective value(s).
 */
struct transform_obj_penalize_t final : public coco_problem_transformed_t {
  const double factor;

  transform_obj_penalize_t(coco_problem_t *inner_problem, const double factor)
      : coco_problem_transformed_t(inner_problem, "transform_obj_penalize"), factor(factor) {
    /* No need to update the best value as the best parameter is feasible */
  }

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(const double *x, double *y) override {
    const double *lower_bounds = smallest_values_of_interest.data();
    const double *upper_bounds = largest_values_of_interest.data();
    double penalty = 0.0;
    size_t i;

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }

    for (i = 0; i < number_of_variables; ++i) {
      const double c1 = x[i] - upper_bounds[i];
      const double c2 = lower_bounds[i] - x[i];
      assert(lower_bounds[i] < upper_bounds[i]);
      if (c1 > 0.0) {
        penalty += c1 * c1;
      } else if (c2 > 0.0) {
        penalty += c2 * c2;
      }
    }
    assert(inner_problem != NULL);
    coco_evaluate_function(inner_problem, x, y);

    for (i = 0; i < number_of_objectives; ++i) {
      y[i] += factor * penalty;
    }
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_obj_penalize(coco_problem_t *inner_problem, const double factor) {
  assert(inner_problem != NULL);
  return new transform_obj_penalize_t(inner_problem, factor);
}
