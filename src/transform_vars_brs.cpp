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
class transform_vars_brs_t final : public coco_problem_transformed_t {
public:
  explicit transform_vars_brs_t(coco_problem_ptr_t inner)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_brs"),
        scaled_x(inner_problem->get_number_of_variables()), factors(inner_problem->get_number_of_variables()) {
    /* Function documentation says we should compute 10^(0.5 *
     * (i-1)/(D-1)). Instead we compute the equivalent
     * sqrt(10)^((i-1)/(D-1)) just like the legacy code.
     */
    for (size_t i = 0; i < number_of_variables; ++i)
      factors[i] = pow(sqrt(10.0), static_cast<double>(i) / (static_cast<double>(number_of_variables) - 1.0));
    update_best_parameter();
  }

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    for (size_t i = 0; i < number_of_variables; ++i) {
      double factor = factors[i];
      /* Documentation specifies odd indices and starts indexing
       * from 1, we use all even indices since C starts indexing
       * with 0.
       */
      if (x[i] > 0.0 && i % 2 == 0) {
        factor *= 10.0;
      }
      scaled_x[i] = factor * x[i];
    }
    inner_problem->evaluate_function(scaled_x, y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /** @brief Divides by the factor; the factor is positive, so inner_x has the sign of x. */
  bool inverse(std::span<const double> inner_x, std::span<double> x) const override {
    for (size_t i = 0; i < number_of_variables; ++i) {
      double factor = factors[i];
      if (inner_x[i] > 0.0 && i % 2 == 0) {
        factor *= 10.0;
      }
      x[i] = inner_x[i] / factor;
    }
    return true;
  }

private:
  std::vector<double> scaled_x; /**< @brief The scaled decision values */
  std::vector<double> factors;  /**< @brief sqrt(10)^(i / (n - 1)) */
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_brs(coco_problem_ptr_t inner_problem) {
  return std::make_unique<transform_vars_brs_t>(std::move(inner_problem));
}
