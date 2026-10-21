/**
 * @file transform_vars_enforce_box.cpp
 * @brief Enforce box constraint, returning NaN if box is violated.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"

namespace coco {

/**
 * @brief Returns Inf for solutions outside of the region of interest.
 */
class TransformVarsEnforceBox final : public TransformedProblem {
public:
  explicit TransformVarsEnforceBox(ProblemPtr inner)
      : TransformedProblem(std::move(inner), "transform_vars_enforce_box") {}

  /**
   * @brief Evaluates the transformed objective function.
   */
  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    for (size_t i = 0; i < inner_problem_->number_of_variables(); ++i) {
      if ((x[i] < smallest_values_of_interest_[i]) || (x[i] > largest_values_of_interest_[i])) {
        /* Box constraint is violated, return Inf */
        VectorSetToInf(y);
        return;
      }
    }

    inner_problem_->EvaluateFunction(x, y);

    AssertNotBelowBestValue(x, y);
  }

  /** @brief The decision values are not transformed. */
  bool Inverse(std::span<const double> inner_x, std::span<double> x) const override {
    return IdentityInverse(inner_x, x);
  }
};

} // namespace coco
