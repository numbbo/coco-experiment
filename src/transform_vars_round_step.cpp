/**
 * @file transform_vars_round_step.cpp
 * @brief Implementation of rounding the variables for the step-ellipsoid function
 * TODO: should this be a helper function instead?
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "coco_utilities.cpp"

/**
 * @brief Data that the large-scale step ellipsoid function shares with transform_vars_round_step.
 */
struct f_step_ellipsoid_shared_data_t {
  double zhat_1; /**< @brief contains the value of \hat{z}_1 that is used to compute the fintess */
};

/**
 * @brief Rounds the decision values for the step-ellipsoid function.
 *
 * Also sets \hat{z}_1 in the data shared with the step ellipsoid function.
 */
struct transform_vars_round_step_t final : public coco_problem_transformed_t {
  const double alpha;
  std::vector<double> rounded_x;
  f_step_ellipsoid_shared_data_t *const shared; /**< @brief Owned by the step ellipsoid function */

  transform_vars_round_step_t(coco_problem_ptr_t inner, const double alpha, f_step_ellipsoid_shared_data_t *shared)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_round_step"), alpha(alpha),
        rounded_x(inner_problem->number_of_variables + 1), shared(shared) {
    /* Compute best parameter */
    for (size_t i = 0; i < number_of_variables; i++) {
      if (fabs(best_parameter[i]) > 0.5) {
        best_parameter[i] = coco_double_round(best_parameter[i]);
      } else {
        best_parameter[i] = coco_double_round(alpha * best_parameter[i]) / alpha;
      }
    }
  }

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(const double *x, double *y) override {

    /* multiplication by d to counter-balance the normalization by d*/
    shared->zhat_1 =
        fabs(x[0]) *
        static_cast<double>(
            inner_problem->number_of_variables); /* TODO: Discuss: consider not pre-imptively multiplying by dim to not
                                                   change the outcome of the max in the core function even though we
                                                   might want to keep it as it is since otherwise, the sum part of the
                                                   max may take over as dim increases */
    for (size_t i = 0; i < inner_problem->number_of_variables; ++i) {
      if (fabs(x[i]) > 0.5) {
        rounded_x[i] = coco_double_round(x[i]);
      } else {
        rounded_x[i] = coco_double_round(alpha * x[i]) / alpha;
      }
    }
    inner_problem->evaluate_function(rounded_x.data(), y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /** @brief Rounding is not invertible. */
  bool inverse(const double *, double *) const override { return false; }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_round_step(coco_problem_ptr_t inner_problem, const double alpha,
                                                    f_step_ellipsoid_shared_data_t *shared) {
  return std::make_unique<transform_vars_round_step_t>(std::move(inner_problem), alpha, shared);
}
