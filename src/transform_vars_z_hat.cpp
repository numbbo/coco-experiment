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
class transform_vars_z_hat_t final : public coco_problem_transformed_t {
public:
  transform_vars_z_hat_t(coco_problem_ptr_t inner, std::span<const double> xopt)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_z_hat"), xopt(xopt.begin(), xopt.end()),
        z(inner_problem->get_number_of_variables()) {
    /* Keeps the best parameter of the inner problem, see f_schwefel_bbob_problem_allocate(): at the optimum of
     * the Schwefel functions, the only ones that use this transformation, x[i - 1] = 2 |xopt[i - 1]| and the
     * transformation is the identity. */
  }

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {

    z[0] = x[0];
    for (size_t i = 1; i < number_of_variables; ++i) {
      z[i] = x[i] + 0.25 * (x[i - 1] - 2.0 * fabs(xopt[i - 1]));
    }
    inner_problem->evaluate_function(z, y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /** @brief Solves for x by forward substitution. */
  bool inverse(std::span<const double> inner_x, std::span<double> x) const override {
    if (number_of_variables == 0)
      return true;
    x[0] = inner_x[0];
    for (size_t i = 1; i < number_of_variables; ++i)
      x[i] = inner_x[i] - 0.25 * (x[i - 1] - 2.0 * fabs(xopt[i - 1]));
    return true;
  }

private:
  const std::vector<double> xopt;
  std::vector<double> z; /**< @brief The transformed decision values */
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_z_hat(coco_problem_ptr_t inner_problem, std::span<const double> xopt) {
  return std::make_unique<transform_vars_z_hat_t>(std::move(inner_problem), xopt);
}
