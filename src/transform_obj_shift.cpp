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
class transform_obj_shift_t final : public coco_problem_transformed_t {
public:
  transform_obj_shift_t(coco_problem_ptr_t inner, const double offset)
      : coco_problem_transformed_t(std::move(inner), "transform_obj_shift"), offset(offset) {
    for (size_t i = 0; i < number_of_objectives; i++)
      best_value[i] += offset;
  }

  void evaluate_function(const double *x, double *y) override {
    if (number_of_objectives == 0) {
      inner_problem->evaluate_function(x, y);
      return;
    }

    inner_problem->evaluate_function(x, y);

    for (size_t i = 0; i < number_of_objectives; i++)
      y[i] += offset;

    assert_not_below_best_value(x, y);
  }

  void evaluate_gradient(const double *x, double *y) override {
    inner_problem->evaluate_gradient(x, y);
  }

  /** @brief The decision values are not transformed. */
  bool inverse(const double *inner_x, double *x) const override { return identity_inverse(inner_x, x); }

private:
  const double offset;
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_obj_shift(coco_problem_ptr_t inner_problem, const double offset) {
  return std::make_unique<transform_obj_shift_t>(std::move(inner_problem), offset);
}
