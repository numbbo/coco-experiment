/**
 * @file transform_vars_conditioning.cpp
 * @brief Implementation of conditioning decision values.
 */

#include <math.h>
#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "transform_functions.cpp"

/**
 * @brief Conditions the decision values.
 */
class transform_vars_conditioning_t final : public coco_problem_transformed_t {
public:
  transform_vars_conditioning_t(coco_problem_ptr_t inner, const double alpha, const bool keep_best_parameter)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_conditioning"),
        conditioning(alpha, inner_problem->get_number_of_variables()),
        conditioned_x(inner_problem->get_number_of_variables()) {
    if (!keep_best_parameter)
      update_best_parameter();
  }

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(const double *x, double *y) override {
    conditioning.apply(x, conditioned_x.data());
    inner_problem->evaluate_function(conditioned_x.data(), y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the gradient of the transformed function.
   */
  void evaluate_gradient(const double *x, double *y) override {
    conditioning.apply(x, conditioned_x.data());
    inner_problem->evaluate_gradient(conditioned_x.data(), y);
    conditioning.apply(y, y);
  }

  /** @brief Divides by the scaling coefficients. */
  bool inverse(const double *inner_x, double *x) const override {
    conditioning.inverse(inner_x, x);
    return true;
  }

private:
  const vars_conditioning_t conditioning;
  std::vector<double> conditioned_x; /**< @brief The conditioned decision values */
};

/**
 * @brief Creates the transformation.
 *
 * If keep_best_parameter is set, the best parameter of the inner problem is kept instead of being computed with
 * inverse(). Only the Schwefel functions use this, see f_schwefel_bbob_problem_allocate().
 */
static coco_problem_ptr_t transform_vars_conditioning(coco_problem_ptr_t inner_problem, const double alpha,
                                                      const bool keep_best_parameter = false) {
  return std::make_unique<transform_vars_conditioning_t>(std::move(inner_problem), alpha, keep_best_parameter);
}
