/**
 * @file transform_obj_oscillate.cpp
 * @brief Implementation of oscillating the objective value.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "coco_problem.cpp"

/**
 * @brief Oscillates the objective value(s) of the inner problem.
 */
class transform_obj_oscillate_t final : public coco_problem_transformed_t {
public:
  explicit transform_obj_oscillate_t(coco_problem_ptr_t inner)
      : coco_problem_transformed_t(std::move(inner), "transform_obj_oscillate") {
    /* Compute best value */
    /* Maybe not the most efficient solution */
    evaluate_function(best_parameter.data(), best_value.data());
  }

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(const double *x, double *y) override {
    static const double factor = 0.1;

    inner_problem->evaluate_function(x, y);

    for (size_t i = 0; i < number_of_objectives; i++) {
      if (y[i] != 0) {
        double log_y = log(fabs(y[i])) / factor;
        if (y[i] > 0) {
          y[i] = pow(exp(log_y + 0.49 * (sin(log_y) + sin(0.79 * log_y))), factor);
        } else {
          y[i] = -pow(exp(log_y + 0.49 * (sin(0.55 * log_y) + sin(0.31 * log_y))), factor);
        }
      }
    }
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /** @brief The decision values are not transformed. */
  bool inverse(const double *inner_x, double *x) const override { return identity_inverse(inner_x, x); }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_obj_oscillate(coco_problem_ptr_t inner_problem) {
  return std::make_unique<transform_obj_oscillate_t>(std::move(inner_problem));
}
