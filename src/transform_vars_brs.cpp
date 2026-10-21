/**
 * @file transform_vars_brs.cpp
 * @brief Implementation of the ominous 's_i scaling' of the BBOB Bueche-Rastrigin problem.
 */

#include <math.h>
#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"

/**
 * @brief The 's_i scaling' of the BBOB Bueche-Rastrigin problem.
 */
struct transform_vars_brs_t final : public coco_problem_transformed_t {
  std::vector<double> scaled_x; /**< @brief The scaled decision values */

  explicit transform_vars_brs_t(coco_problem_ptr_t inner)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_brs"),
        scaled_x(inner_problem->number_of_variables) {
    if (coco_problem_best_parameter_not_zero(inner_problem.get())) {
      coco_warning("transform_vars_brs(): 'best_parameter' not updated, set to NAN");
      coco_vector_set_to_nan(inner_problem->best_parameter.data(), inner_problem->number_of_variables);
    }
  }

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(const double *x, double *y) override {
    size_t i;
    double factor;

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }

    for (i = 0; i < number_of_variables; ++i) {
      /* Function documentation says we should compute 10^(0.5 *
       * (i-1)/(D-1)). Instead we compute the equivalent
       * sqrt(10)^((i-1)/(D-1)) just like the legacy code.
       */
      factor = pow(sqrt(10.0), (double)(long)i / ((double)(long)number_of_variables - 1.0));
      /* Documentation specifies odd indices and starts indexing
       * from 1, we use all even indices since C starts indexing
       * with 0.
       */
      if (x[i] > 0.0 && i % 2 == 0) {
        factor *= 10.0;
      }
      scaled_x[i] = factor * x[i];
    }
    coco_evaluate_function(inner_problem.get(), scaled_x.data(), y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_brs(coco_problem_ptr_t inner_problem) {
  return std::make_unique<transform_vars_brs_t>(std::move(inner_problem));
}
