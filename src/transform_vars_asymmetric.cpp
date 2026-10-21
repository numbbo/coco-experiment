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
#include "brentq.cpp"

/**
 * @brief Data type for univariate function tasy_uv.
 */
struct tasy_data {
  double beta;
  size_t i;
  size_t n;
};

/**
 * @brief Univariate asymmetric non-linear transformation.
 */
static double tasy_uv(double xi, const tasy_data *d) {
  double yi;
  double exponent;
  if (xi > 0.0) {
    exponent = 1.0 + ((d->beta * static_cast<double>(d->i)) / (static_cast<double>(d->n) - 1.0)) * sqrt(xi);
    yi = pow(xi, exponent);
  } else {
    yi = xi;
  }
  return yi;
}

/**
 * @brief Inverse of asymmetric non-linear transformation tasy_uv obtained with brentq.
 */
static double tasy_uv_inv(double yi, const tasy_data *d) {
  return brentinv([d](double xi) { return tasy_uv(xi, d); }, yi);
}

/**
 * @brief Applies the asymmetric non-linear transformation to the decision values.
 */
class transform_vars_asymmetric_t final : public coco_problem_transformed_t {
public:
  transform_vars_asymmetric_t(coco_problem_ptr_t inner, const double beta)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_asymmetric"),
        asymmetric_x(inner_problem->get_number_of_variables()), beta(beta) {
    if (inner_problem->get_number_of_objectives() > 0 && inner_problem->best_parameter_not_zero()) {
      coco_warning("transform_vars_asymmetric(): 'best_parameter' not updated, set to NAN");
      inner_problem->set_best_parameter_to_nan();
    }
  }

  /**
   * @brief Multivariate, coordinate-wise, asymmetric non-linear transformation.
   */
  void tasy(const double *x) {
    tasy_data d;

    d.beta = beta;
    d.n = number_of_variables;

    for (size_t i = 0; i < number_of_variables; ++i) {
      d.i = i;
      asymmetric_x[i] = tasy_uv(x[i], &d);
    }
  }

  /**
   * @brief Evaluates the transformed function.
   */
  void evaluate_function(const double *x, double *y) override {
    if (inner_problem->get_number_of_objectives() == 0) {
      inner_problem->evaluate_function(x, y);
      return;
    }

    tasy(x);

    inner_problem->evaluate_function(asymmetric_x.data(), y);

    assert_not_below_best_value(x, y);
  }

  /**
   * @brief Evaluates the transformed constraint.
   */
  void evaluate_constraint(const double *x, double *y, int update_counter) override {
    double exponent;

    if (inner_problem->get_number_of_constraints() == 0) {
      inner_problem->evaluate_constraint(x, y, update_counter);
      return;
    }

    /* FIXME (OME): Old  pre-logger
    tasy(x);

    coco_evaluate_constraint(inner_problem.get(), asymmetric_x, y);
    */
    for (size_t i = 0; i < number_of_variables; ++i) {
      if (x[i] > 0.0) {
        exponent =
            1.0 + ((beta * static_cast<double>(i)) / (static_cast<double>(number_of_variables) - 1.0)) * sqrt(x[i]);
        asymmetric_x[i] = pow(x[i], exponent);
      } else {
        asymmetric_x[i] = x[i];
      }
    }
    inner_problem->evaluate_constraint(asymmetric_x.data(), y, update_counter);
  }

  /** @brief Inverts tasy numerically (with Brent's method). */
  bool inverse(const double *inner_x, double *x) const override {
    tasy_data d;
    d.beta = beta;
    d.n = number_of_variables;
    for (size_t i = 0; i < number_of_variables; ++i) {
      d.i = i;
      x[i] = tasy_uv_inv(inner_x[i], &d);
    }
    return true;
  }

private:
  std::vector<double> asymmetric_x; /**< @brief The transformed decision values */
  const double beta;
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_asymmetric(coco_problem_ptr_t inner_problem, const double beta) {
  return std::make_unique<transform_vars_asymmetric_t>(std::move(inner_problem), beta);
}
