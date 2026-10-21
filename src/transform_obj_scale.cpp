/**
 * @file transform_obj_scale.cpp
 * @brief Implementation of scaling the objective value by the given factor.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"

namespace coco {

/**
 * @brief Scales the objective value(s) of the inner problem by a factor.
 */
class TransformObjScale final : public TransformedProblem {
public:
  TransformObjScale(ProblemPtr inner, const double factor)
      : TransformedProblem(std::move(inner), "transform_obj_scale"), factor_(factor) {
    for (size_t i = 0; i < number_of_objectives_; ++i)
      best_value_[i] *= factor;
  }

  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    if (number_of_objectives_ == 0) {
      inner_problem_->EvaluateFunction(x, y);
      return;
    }

    inner_problem_->EvaluateFunction(x, y);

    for (size_t i = 0; i < number_of_objectives_; i++)
      y[i] *= factor_;

    AssertNotBelowBestValue(x, y);
  }

  void EvaluateGradient(std::span<const double> x, std::span<double> y) override {

    inner_problem_->EvaluateGradient(x, y);

    for (size_t i = 0; i < number_of_variables_; ++i) {
      y[i] *= factor_;
    }
  }

  /** @brief The decision values are not transformed. */
  bool Inverse(std::span<const double> inner_x, std::span<double> x) const override {
    return IdentityInverse(inner_x, x);
  }

private:
  const double factor_;
};

} // namespace coco
