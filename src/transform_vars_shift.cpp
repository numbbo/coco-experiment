/**
 * @file transform_vars_shift.cpp
 * @brief Implementation of shifting all decision values by an offset.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "transform_functions.cpp"

namespace coco {

/**
 * @brief Shifts all decision values by an offset.
 */
class TransformVarsShift final : public TransformedProblem {
public:
  /**
   * CAVEAT: when shifting the constraint only, the best_value of best_parameter
   *         will get in an inconsistent state.
   */
  TransformVarsShift(ProblemPtr inner, std::span<const double> offset, const bool shift_constraint_only)
      : TransformedProblem(std::move(inner), "transform_vars_shift"), shift_(offset),
        shifted_x_(inner_problem_->number_of_variables()), shift_constraint_only_(shift_constraint_only) {
    assert(offset.size() == number_of_variables_);

    UpdateBestParameter();

    /* Update the initial solution if any */
    if (!initial_solution_.empty())
      shift_.Inverse(initial_solution_, initial_solution_);
  }

  /**
   * @brief Evaluates the transformed objective function.
   */
  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    if (inner_problem_->number_of_objectives() == 0 || shift_constraint_only_ != 0) {
      inner_problem_->EvaluateFunction(x, y);
      return;
    }

    shift_.Apply(x, shifted_x_);
    inner_problem_->EvaluateFunction(shifted_x_, y);

    AssertNotBelowBestValue(x, y);
  }

  /**
   * @brief Evaluates the transformed constraint function.
   */
  void EvaluateConstraint(std::span<const double> x, std::span<double> y, int update_counter) override {

    if (inner_problem_->number_of_constraints() == 0) {
      inner_problem_->EvaluateConstraint(x, y, update_counter);
      return;
    }

    shift_.Apply(x, shifted_x_);
    inner_problem_->EvaluateConstraint(shifted_x_, y, update_counter);
  }

  /**
   * @brief Evaluates the gradient of the transformed function at x
   */
  void EvaluateGradient(std::span<const double> x, std::span<double> y) override {
    shift_.Apply(x, shifted_x_);
    inner_problem_->EvaluateGradient(shifted_x_, y);
  }

  /** @brief Adds the offset again. */
  bool Inverse(std::span<const double> inner_x, std::span<double> x) const override {
    shift_.Inverse(inner_x, x);
    return true;
  }

private:
  const VarsShift shift_;
  std::vector<double> shifted_x_;
  const bool shift_constraint_only_;
};

} // namespace coco
