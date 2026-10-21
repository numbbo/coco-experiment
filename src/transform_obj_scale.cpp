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
      y[i] *= factor;

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

    for (size_t i = 0; i < number_of_variables; ++i) {
      y[i] *= factor;
    }
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_obj_scale(coco_problem_ptr_t inner_problem, const double factor) {
  return std::make_unique<transform_obj_scale_t>(std::move(inner_problem), factor);
}
