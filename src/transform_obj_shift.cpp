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

  transform_obj_shift_t(coco_problem_ptr_t inner, const double offset)
      : coco_problem_transformed_t(std::move(inner), "transform_obj_shift"), offset(offset) {
    for (size_t i = 0; i < number_of_objectives; i++)
      best_value[i] += offset;
  }

  void evaluate_function(const double *x, double *y) override {
    int feasible;

    if (number_of_objectives == 0) {
      coco_evaluate_function(inner_problem.get(), x, y);
      return;
    }

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }

    coco_evaluate_function(inner_problem.get(), x, y);

    for (size_t i = 0; i < number_of_objectives; i++)
      y[i] += offset;

    if (number_of_constraints > 0) {
      std::vector<double> cons_values(number_of_constraints);
      feasible = is_feasible(x, cons_values.data());
      if (feasible)
        assert(y[0] + 1e-13 >= best_value[0]);
    } else
      assert(y[0] + 1e-13 >= best_value[0]);
  }

  void evaluate_gradient(const double *x, double *y) override {
    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }

    inner_problem->evaluate_gradient(x, y);
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_obj_shift(coco_problem_ptr_t inner_problem, const double offset) {
  return std::make_unique<transform_obj_shift_t>(std::move(inner_problem), offset);
}
