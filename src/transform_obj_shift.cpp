/**
 * @file transform_obj_shift.cpp
 * @brief Implementation of shifting the objective value by the given offset.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"

namespace coco {

/**
 * @brief Shifts the objective value(s) of the inner problem by an offset.
 */
class TransformObjShift final : public TransformedProblem {
public:
  TransformObjShift(ProblemPtr inner, const double offset)
      : TransformedProblem(std::move(inner), "transform_obj_shift"), offset_(offset) {
    for (size_t i = 0; i < number_of_objectives_; i++)
      best_value_[i] += offset;
  }

  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    if (number_of_objectives_ == 0) {
      inner_problem_->EvaluateFunction(x, y);
      return;
    }

    inner_problem_->EvaluateFunction(x, y);

    for (size_t i = 0; i < number_of_objectives_; i++)
      y[i] += offset_;

    AssertNotBelowBestValue(x, y);
  }

  void EvaluateGradient(std::span<const double> x, std::span<double> y) override {
    inner_problem_->EvaluateGradient(x, y);
  }

  /** @brief The decision values are not transformed. */
  bool Inverse(std::span<const double> inner_x, std::span<double> x) const override {
    return IdentityInverse(inner_x, x);
  }

private:
  const double offset_;
};

} // namespace coco
