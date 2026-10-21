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

/**
 * @brief Applies the asymmetric non-linear transformation to the decision values.
 */
class transform_vars_asymmetric_t final : public coco_problem_transformed_t {
public:
  transform_vars_asymmetric_t(coco_problem_ptr_t inner, const double beta)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_asymmetric"),
        asymmetric(beta, inner_problem->get_number_of_variables()),
        asymmetric_x(inner_problem->get_number_of_variables()) {
    update_best_parameter();
  }

  /**
   * @brief Evaluates the transformed function.
   */
  void evaluate_function(const double *x, double *y) override {
    if (inner_problem->get_number_of_objectives() == 0) {
      inner_problem->evaluate_function(x, y);
      return;
    }

    asymmetric.apply(x, asymmetric_x.data());
    inner_problem->evaluate_function(asymmetric_x.data(), y);

    assert_not_below_best_value(x, y);
  }

  /**
   * @brief Evaluates the transformed constraint.
   */
  void evaluate_constraint(const double *x, double *y, int update_counter) override {
    if (inner_problem->get_number_of_constraints() == 0) {
      inner_problem->evaluate_constraint(x, y, update_counter);
      return;
    }

    asymmetric.apply(x, asymmetric_x.data());
    inner_problem->evaluate_constraint(asymmetric_x.data(), y, update_counter);
  }

  /** @brief Inverts the transformation numerically. */
  bool inverse(const double *inner_x, double *x) const override {
    asymmetric.inverse(inner_x, x);
    return true;
  }

private:
  const vars_asymmetric_t asymmetric;
  std::vector<double> asymmetric_x; /**< @brief The transformed decision values */
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_asymmetric(coco_problem_ptr_t inner_problem, const double beta) {
  return std::make_unique<transform_vars_asymmetric_t>(std::move(inner_problem), beta);
}
