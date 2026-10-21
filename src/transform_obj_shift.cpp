/**
 * @file transform_obj_shift.cpp
 * @brief Implementation of shifting the objective value by the given offset.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"

/**
 * @brief Shifts the objective value(s) of the inner problem by an offset.
 */
struct transform_obj_shift_t final : public coco_problem_transformed_t {
  const double offset;

  transform_obj_shift_t(coco_problem_t *inner_problem, const double offset)
      : coco_problem_transformed_t(inner_problem, "transform_obj_shift"), offset(offset) {
    size_t i;
    for (i = 0; i < number_of_objectives; i++)
      best_value[i] += offset;
  }

  void evaluate_function(const double *x, double *y) override {
    int is_feasible;
    size_t i;

    if (number_of_objectives == 0) {
      coco_evaluate_function(inner_problem, x, y);
      return;
    }

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }

    coco_evaluate_function(inner_problem, x, y);

    for (i = 0; i < number_of_objectives; i++)
      y[i] += offset;

    if (number_of_constraints > 0) {
      std::vector<double> cons_values(number_of_constraints);
      is_feasible = coco_is_feasible(this, x, cons_values.data());
      if (is_feasible)
        assert(y[0] + 1e-13 >= best_value[0]);
    } else
      assert(y[0] + 1e-13 >= best_value[0]);
  }

  void evaluate_gradient(const double *x, double *y) override {
    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }

    bbob_evaluate_gradient(inner_problem, x, y);
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_obj_shift(coco_problem_t *inner_problem, const double offset) {
  return new transform_obj_shift_t(inner_problem, offset);
}
