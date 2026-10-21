/**
 * @file transform_vars_enforce_box.cpp
 * @brief Enforce box constraint, returning NaN if box is violated.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"

/**
 * @brief Returns Inf for solutions outside of the region of interest.
 */
struct transform_vars_enforce_box_t final : public coco_problem_transformed_t {
  explicit transform_vars_enforce_box_t(coco_problem_ptr_t inner)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_enforce_box") {}

  /**
   * @brief Evaluates the transformed objective function.
   */
  void evaluate_function(const double *x, double *y) override {
    int feasible;

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }

    for (size_t i = 0; i < inner_problem->number_of_variables; ++i) {
      if ((x[i] < smallest_values_of_interest[i]) || (x[i] > largest_values_of_interest[i])) {
        /* Box constraint is violated, return Inf */
        coco_vector_set_to_inf(y, coco_problem_get_number_of_objectives(this));
        return;
      }
    }

    coco_evaluate_function(inner_problem.get(), x, y);

    if (number_of_constraints > 0) {
      std::vector<double> cons_values(number_of_constraints);
      feasible = is_feasible(x, cons_values.data());
      if (feasible) {
        assert(y[0] + 1e-13 >= best_value[0]);
      }
    } else {
      assert(y[0] + 1e-13 >= best_value[0]);
    }
  }

  /** @brief The decision values are not transformed. */
  bool inverse(const double *inner_x, double *x) const override { return identity_inverse(inner_x, x); }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_enforce_box(coco_problem_ptr_t inner_problem) {
  return std::make_unique<transform_vars_enforce_box_t>(std::move(inner_problem));
}
