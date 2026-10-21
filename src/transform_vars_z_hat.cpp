/**
 * @file transform_vars_z_hat.cpp
 * @brief Implementation of the z^hat transformation of decision values for the BBOB Schwefel problem.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"

/**
 * @brief The z^hat transformation of decision values for the BBOB Schwefel problem.
 */
struct transform_vars_z_hat_t final : public coco_problem_transformed_t {
  const std::vector<double> xopt;
  std::vector<double> z; /**< @brief The transformed decision values */

  transform_vars_z_hat_t(coco_problem_ptr_t inner, const double *xopt)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_z_hat"),
        xopt(xopt, xopt + inner_problem->number_of_variables), z(inner_problem->number_of_variables) {
    /* TODO: implement best_parameter transformation if needed in the case of not zero:
       see also issue #814.
    The correct update of best_parameter seems not too difficult and should not anymore
    break the current implementation of the Schwefel function.
     coco_warning("transform_vars_z_hat(): 'best_parameter' not updated");

    This:
    size_t i;
    if (problem->best_parameter != NULL)
          for (i = 1; i < problem->number_of_variables; ++i)
            problem->best_parameter[i] -= 0.25 * (problem->best_parameter[i - 1] - 2.0 * fabs(data->xopt[i - 1]));
    should do, but gives
    COCO INFO: ..., d=2, running: f18.Assertion failed: (about_equal_value(hypervolume, 8.1699208579037619e-05)), function
    test_coco_archive_extreme_solutions, file ./test_coco_archive.cpp, line 123.
    */
    if (strstr(coco_problem_get_id(inner_problem.get()), "schwefel") == nullptr) {
      coco_warning("transform_vars_z_hat(): 'best_parameter' not updated, set to NAN.");
      coco_vector_set_to_nan(best_parameter.data(), number_of_variables);
    }
  }

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(const double *x, double *y) override {
    size_t i;

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }

    z[0] = x[0];
    for (i = 1; i < number_of_variables; ++i) {
      z[i] = x[i] + 0.25 * (x[i - 1] - 2.0 * fabs(xopt[i - 1]));
    }
    coco_evaluate_function(inner_problem.get(), z.data(), y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_z_hat(coco_problem_ptr_t inner_problem, const double *xopt) {
  return std::make_unique<transform_vars_z_hat_t>(std::move(inner_problem), xopt);
}
