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
class transform_vars_enforce_box_t final : public coco_problem_transformed_t {
public:
  explicit transform_vars_enforce_box_t(coco_problem_ptr_t inner)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_enforce_box") {}

  /**
   * @brief Evaluates the transformed objective function.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    for (size_t i = 0; i < inner_problem->get_number_of_variables(); ++i) {
      if ((x[i] < smallest_values_of_interest[i]) || (x[i] > largest_values_of_interest[i])) {
        /* Box constraint is violated, return Inf */
        coco_vector_set_to_inf(y);
        return;
      }
    }

    inner_problem->evaluate_function(x, y);

    assert_not_below_best_value(x, y);
  }

  /** @brief The decision values are not transformed. */
  bool inverse(std::span<const double> inner_x, std::span<double> x) const override {
    return identity_inverse(inner_x, x);
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_enforce_box(coco_problem_ptr_t inner_problem) {
  return std::make_unique<transform_vars_enforce_box_t>(std::move(inner_problem));
}
