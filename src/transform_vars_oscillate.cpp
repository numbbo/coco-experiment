/**
 * @file transform_vars_oscillate.cpp
 * @brief Implementation of oscillating the decision values.
 * @author ??
 * @author Paul Dufosse
 * @note Edited to fulfill needs from the constrained test bed.
 */

#include <math.h>
#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "brentq.cpp"

/**
 * @brief Data type for univariate function tosz_uv
 */
struct tosz_data {
  double alpha;
};

/**
 * @brief Univariate oscillating non-linear transformation.
 */
static double tosz_uv(double xi, const tosz_data *d) {
  double yi;
  double tmp, base;
  if (xi > 0.0) {
    tmp = log(xi) / d->alpha;
    base = exp(tmp + 0.49 * (sin(tmp) + sin(0.79 * tmp)));
    yi = pow(base, d->alpha);
  } else if (xi < 0.0) {
    tmp = log(-xi) / d->alpha;
    base = exp(tmp + 0.49 * (sin(0.55 * tmp) + sin(0.31 * tmp)));
    yi = -pow(base, d->alpha);
  } else {
    yi = 0.0;
  }
  return yi;
}

/**
 * @brief Inverse of oscillating non-linear transformation tosz_uv obtained with brentq.
 */
static double tosz_uv_inv(double yi, const tosz_data *d) {
  return brentinv([d](double xi) { return tosz_uv(xi, d); }, yi);
}

/**
 * @brief Applies the oscillating non-linear transformation to the decision values.
 */
class transform_vars_oscillate_t final : public coco_problem_transformed_t {
public:
  explicit transform_vars_oscillate_t(coco_problem_ptr_t inner)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_oscillate"), alpha(0.1),
        oscillated_x(inner_problem->get_number_of_variables()) {}

  /**
   * @brief Multivariate, coordinate-wise, oscillating non-linear transformation.
   */
  void tosz(const double *x) {
    tosz_data d;

    d.alpha = alpha;

    for (size_t i = 0; i < number_of_variables; ++i) {
      oscillated_x[i] = tosz_uv(x[i], &d);
    }
  }

  /**
   * @brief Evaluates the transformed objective functions.
   */
  void evaluate_function(const double *x, double *y) override {
    if (inner_problem->get_number_of_objectives() == 0) {
      inner_problem->evaluate_function(x, y);
      return;
    }

    tosz(x);

    inner_problem->evaluate_function(oscillated_x.data(), y);

    assert_not_below_best_value(x, y);
  }

  /**
   * @brief Evaluates the transformed constraints.
   */
  void evaluate_constraint(const double *x, double *y, int update_counter) override {
    double tmp, base;

    if (inner_problem->get_number_of_constraints() == 0) {
      inner_problem->evaluate_constraint(x, y, update_counter);
      return;
    }

    for (size_t i = 0; i < number_of_variables; ++i) {
      if (x[i] > 0.0) {
        tmp = log(x[i]) / alpha;
        base = exp(tmp + 0.49 * (sin(tmp) + sin(0.79 * tmp)));
        oscillated_x[i] = pow(base, alpha);
      } else if (x[i] < 0.0) {
        tmp = log(-x[i]) / alpha;
        base = exp(tmp + 0.49 * (sin(0.55 * tmp) + sin(0.31 * tmp)));
        oscillated_x[i] = -pow(base, alpha);
      } else {
        oscillated_x[i] = 0.0;
      }
    }
    inner_problem->evaluate_constraint(oscillated_x.data(), y, update_counter);
  }

  /** @brief Inverts tosz numerically (with Brent's method). */
  bool inverse(const double *inner_x, double *x) const override {
    tosz_data d;
    d.alpha = alpha;
    for (size_t i = 0; i < number_of_variables; ++i)
      x[i] = tosz_uv_inv(inner_x[i], &d);
    return true;
  }

private:
  const double alpha;
  std::vector<double> oscillated_x;
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_oscillate(coco_problem_ptr_t inner_problem) {
  return std::make_unique<transform_vars_oscillate_t>(std::move(inner_problem));
}
