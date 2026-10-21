/**
 * @file transform_vars_asymmetric.cpp
 * @brief Implementation of performing an asymmetric transformation on decision values.
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
 * @brief Applies the asymmetric non-linear transformation to the decision values.
 */
class TransformVarsAsymmetric final : public TransformedProblem {
public:
  TransformVarsAsymmetric(ProblemPtr inner, const double beta)
      : TransformedProblem(std::move(inner), "transform_vars_asymmetric"),
        asymmetric_(beta, inner_problem_->number_of_variables()), asymmetric_x_(inner_problem_->number_of_variables()) {
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

    asymmetric_.Apply(x, asymmetric_x_);
    inner_problem_->EvaluateFunction(asymmetric_x_, y);

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

    asymmetric_.Apply(x, asymmetric_x_);
    inner_problem_->EvaluateConstraint(asymmetric_x_, y, update_counter);
  }

  /** @brief Inverts the transformation numerically. */
  bool Inverse(std::span<const double> inner_x, std::span<double> x) const override {
    asymmetric_.Inverse(inner_x, x);
    return true;
  }

private:
  const VarsAsymmetric asymmetric_;
  std::vector<double> asymmetric_x_; /**< @brief The transformed decision values */
};

} // namespace coco
