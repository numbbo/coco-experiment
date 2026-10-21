/**
 * @file transform_obj_penalize.cpp
 * @brief Implementation of adding a penalty to the objective value for solutions outside of the ROI in the
 * decision space.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "transform_functions.cpp"

/**
 * @brief Adds a penalty for solutions outside of the region of interest to the objective value(s).
 */
class transform_obj_penalize_t final : public coco_problem_transformed_t {
public:
  transform_obj_penalize_t(coco_problem_ptr_t inner, const double factor)
      : coco_problem_transformed_t(std::move(inner), "transform_obj_penalize"), factor(factor) {
    /* No need to update the best value as the best parameter is feasible */
  }

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    const double penalty = coco_boundary_penalty(x, smallest_values_of_interest, largest_values_of_interest);
    inner_problem->evaluate_function(x, y);

    for (size_t i = 0; i < number_of_objectives; ++i) {
      y[i] += factor * penalty;
    }
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /** @brief The decision values are not transformed. */
  bool inverse(std::span<const double> inner_x, std::span<double> x) const override {
    return identity_inverse(inner_x, x);
  }

private:
  const double factor;
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_obj_penalize(coco_problem_ptr_t inner_problem, const double factor) {
  return std::make_unique<transform_obj_penalize_t>(std::move(inner_problem), factor);
}
