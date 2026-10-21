/**
 * @file transform_vars_shift.cpp
 * @brief Implementation of shifting all decision values by an offset.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "transform_functions.cpp"

/**
 * @brief Shifts all decision values by an offset.
 */
class transform_vars_shift_t final : public coco_problem_transformed_t {
public:
  /**
   * CAVEAT: when shifting the constraint only, the best_value of best_parameter
   *         will get in an inconsistent state.
   */
  transform_vars_shift_t(coco_problem_ptr_t inner, std::span<const double> offset, const bool shift_constraint_only)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_shift"), shift(offset),
        shifted_x(inner_problem->get_number_of_variables()), shift_constraint_only(shift_constraint_only) {
    assert(offset.size() == number_of_variables);

    update_best_parameter();

    /* Update the initial solution if any */
    if (!initial_solution.empty())
      shift.inverse(initial_solution, initial_solution);
  }

  /**
   * @brief Evaluates the transformed objective function.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    if (inner_problem->get_number_of_objectives() == 0 || shift_constraint_only != 0) {
      inner_problem->evaluate_function(x, y);
      return;
    }

    shift.apply(x, shifted_x);
    inner_problem->evaluate_function(shifted_x, y);

    assert_not_below_best_value(x, y);
  }

  /**
   * @brief Evaluates the transformed constraint function.
   */
  void evaluate_constraint(std::span<const double> x, std::span<double> y, int update_counter) override {

    if (inner_problem->get_number_of_constraints() == 0) {
      inner_problem->evaluate_constraint(x, y, update_counter);
      return;
    }

    shift.apply(x, shifted_x);
    inner_problem->evaluate_constraint(shifted_x, y, update_counter);
  }

  /**
   * @brief Evaluates the gradient of the transformed function at x
   */
  void evaluate_gradient(std::span<const double> x, std::span<double> y) override {
    shift.apply(x, shifted_x);
    inner_problem->evaluate_gradient(shifted_x, y);
  }

  /** @brief Adds the offset again. */
  bool inverse(std::span<const double> inner_x, std::span<double> x) const override {
    shift.inverse(inner_x, x);
    return true;
  }

private:
  const vars_shift_t shift;
  std::vector<double> shifted_x;
  const bool shift_constraint_only;
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_shift(coco_problem_ptr_t inner_problem, std::span<const double> offset,
                                               const bool shift_constraint_only) {
  return std::make_unique<transform_vars_shift_t>(std::move(inner_problem), offset, shift_constraint_only);
}
