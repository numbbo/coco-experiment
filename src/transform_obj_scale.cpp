/**
 * @file transform_obj_scale.cpp
 * @brief Implementation of scaling the objective value by the given factor.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"

/**
 * @brief Scales the objective value(s) of the inner problem by a factor.
 */
struct transform_obj_scale_t final : public coco_problem_transformed_t {
  const double factor;

  transform_obj_scale_t(coco_problem_ptr_t inner, const double factor)
      : coco_problem_transformed_t(std::move(inner), "transform_obj_scale"), factor(factor) {
    for (size_t i = 0; i < number_of_objectives; ++i)
      best_value[i] *= factor;
  }

  void evaluate_function(const double *x, double *y) override {
    if (number_of_objectives == 0) {
      inner_problem->evaluate_function(x, y);
      return;
    }

    inner_problem->evaluate_function(x, y);

    for (size_t i = 0; i < number_of_objectives; i++)
      y[i] *= factor;

    assert_not_below_best_value(x, y);
  }

  void evaluate_gradient(const double *x, double *y) override {

    inner_problem->evaluate_gradient(x, y);

    for (size_t i = 0; i < number_of_variables; ++i) {
      y[i] *= factor;
    }
  }

  /** @brief The decision values are not transformed. */
  bool inverse(const double *inner_x, double *x) const override { return identity_inverse(inner_x, x); }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_obj_scale(coco_problem_ptr_t inner_problem, const double factor) {
  return std::make_unique<transform_obj_scale_t>(std::move(inner_problem), factor);
}
