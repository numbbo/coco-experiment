/**
 * @file transform_obj_power.cpp
 * @brief Implementation of raising the objective value to the power of a given exponent.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "coco_problem.cpp"

/**
 * @brief Raises the objective value(s) of the inner problem to the power of an exponent.
 */
class transform_obj_power_t final : public coco_problem_transformed_t {
public:
  transform_obj_power_t(coco_problem_ptr_t inner, const double exponent)
      : coco_problem_transformed_t(std::move(inner), "transform_obj_power"), exponent(exponent) {
    /* Compute best value */
    evaluate_function(best_parameter, best_value);
  }

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {

    inner_problem->evaluate_function(x, y);

    for (size_t i = 0; i < number_of_objectives; i++) {
      y[i] = pow(y[i], exponent);
    }
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /** @brief The decision values are not transformed. */
  bool inverse(std::span<const double> inner_x, std::span<double> x) const override {
    return identity_inverse(inner_x, x);
  }

private:
  const double exponent;
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_obj_power(coco_problem_ptr_t inner_problem, const double exponent) {
  return std::make_unique<transform_obj_power_t>(std::move(inner_problem), exponent);
}
