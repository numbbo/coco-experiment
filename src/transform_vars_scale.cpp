/**
 * @file transform_vars_scale.cpp
 * @brief Implementation of scaling decision values by a given factor.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"

/**
 * @brief Scales the decision values by a factor.
 */
struct transform_vars_scale_t final : public coco_problem_transformed_t {
  const double factor;
  std::vector<double> scaled_x; /**< @brief The scaled decision values */

  transform_vars_scale_t(coco_problem_ptr_t inner, const double factor)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_scale"), factor(factor),
        scaled_x(inner_problem->number_of_variables) {
    /* Compute best parameter */
    if (factor != 0.) {
      for (size_t i = 0; i < number_of_variables; i++) {
        best_parameter[i] /= factor;
      }
    } /* else error? */
  }

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(const double *x, double *y) override {

    for (size_t i = 0; i < number_of_variables; ++i) {
      scaled_x[i] = factor * x[i];
    }
    inner_problem->evaluate_function(scaled_x.data(), y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /** @brief Divides by the factor. */
  bool inverse(const double *inner_x, double *x) const override {
    if (factor == 0.)
      return false;
    for (size_t i = 0; i < number_of_variables; ++i)
      x[i] = inner_x[i] / factor;
    return true;
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_scale(coco_problem_ptr_t inner_problem, const double factor) {
  return std::make_unique<transform_vars_scale_t>(std::move(inner_problem), factor);
}
