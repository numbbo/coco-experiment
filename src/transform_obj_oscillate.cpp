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
struct transform_obj_oscillate_t final : public coco_problem_transformed_t {
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

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }
    coco_evaluate_function(inner_problem.get(), x, y);

    for (size_t i = 0; i < number_of_objectives; i++) {
      if (y[i] != 0) {
        double log_y;
        log_y = log(fabs(y[i])) / factor;
        if (y[i] > 0) {
          y[i] = pow(exp(log_y + 0.49 * (sin(log_y) + sin(0.79 * log_y))), factor);
        } else {
          y[i] = -pow(exp(log_y + 0.49 * (sin(0.55 * log_y) + sin(0.31 * log_y))), factor);
        }
      }
    }
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_obj_oscillate(coco_problem_ptr_t inner_problem) {
  return std::make_unique<transform_obj_oscillate_t>(std::move(inner_problem));
}
