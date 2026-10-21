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

namespace coco {

/**
 * @brief Applies the oscillating non-linear transformation to the decision values.
 */
class TransformVarsOscillate final : public TransformedProblem {
public:
  explicit TransformVarsOscillate(ProblemPtr inner)
      : TransformedProblem(std::move(inner), "transform_vars_oscillate"),
        oscillated_x_(inner_problem_->number_of_variables()) {
    UpdateBestParameter();
  }

  /**
   * @brief Evaluates the transformed function.
   */
  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    if (inner_problem_->number_of_objectives() == 0) {
      inner_problem_->EvaluateFunction(x, y);
      return;
    }

    oscillate_.Apply(x, oscillated_x_);
    inner_problem_->EvaluateFunction(oscillated_x_, y);

    AssertNotBelowBestValue(x, y);
  }

  /**
   * @brief Evaluates the transformed constraint.
   */
  void EvaluateConstraint(std::span<const double> x, std::span<double> y, int update_counter) override {
    if (inner_problem_->number_of_constraints() == 0) {
      inner_problem_->EvaluateConstraint(x, y, update_counter);
      return;
    }

    oscillate_.Apply(x, oscillated_x_);
    inner_problem_->EvaluateConstraint(oscillated_x_, y, update_counter);
  }

  /** @brief Inverts the oscillation numerically. */
  bool Inverse(std::span<const double> inner_x, std::span<double> x) const override {
    oscillate_.Inverse(inner_x, x);
    return true;
  }

private:
  const VarsOscillate oscillate_{};
  std::vector<double> oscillated_x_;
};

} // namespace coco
