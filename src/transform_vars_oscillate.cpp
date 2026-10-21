/**
 * @file transform_vars_oscillate.cpp
 * @brief Implementation of oscillating the decision values.
 * @author ??
 * @author Paul Dufosse
 * @note Edited to fulfill needs from the constrained test bed.
 */

#include <math.h>
#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "transform_functions.cpp"

/**
 * @brief Applies the oscillating non-linear transformation to the decision values.
 */
class transform_vars_oscillate_t final : public coco_problem_transformed_t {
public:
  explicit transform_vars_oscillate_t(coco_problem_ptr_t inner)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_oscillate"),
        oscillated_x(inner_problem->get_number_of_variables()) {
    update_best_parameter();
  }

  /**
   * @brief Evaluates the transformed function.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    if (inner_problem->get_number_of_objectives() == 0) {
      inner_problem->evaluate_function(x, y);
      return;
    }

    oscillate.apply(x, oscillated_x);
    inner_problem->evaluate_function(oscillated_x, y);

    assert_not_below_best_value(x, y);
  }

  /**
   * @brief Evaluates the transformed constraint.
   */
  void evaluate_constraint(std::span<const double> x, std::span<double> y, int update_counter) override {
    if (inner_problem->get_number_of_constraints() == 0) {
      inner_problem->evaluate_constraint(x, y, update_counter);
      return;
    }

    oscillate.apply(x, oscillated_x);
    inner_problem->evaluate_constraint(oscillated_x, y, update_counter);
  }

  /** @brief Inverts the oscillation numerically. */
  bool inverse(std::span<const double> inner_x, std::span<double> x) const override {
    oscillate.inverse(inner_x, x);
    return true;
  }

private:
  const vars_oscillate_t oscillate{};
  std::vector<double> oscillated_x;
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_oscillate(coco_problem_ptr_t inner_problem) {
  return std::make_unique<transform_vars_oscillate_t>(std::move(inner_problem));
}
