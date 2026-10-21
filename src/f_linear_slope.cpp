/**
 * @file f_linear_slope.cpp
 * @brief Implementation of the linear slope function and problem.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_norm_by_dim.cpp"
#include "transform_obj_shift.cpp"

/**
 * @brief Implements the linear slope function without connections to any COCO
 * structures.
 */
static double f_linear_slope_raw(const double *x, const size_t number_of_variables, const double *best_parameter) {

  static const double alpha = 100.0;
  size_t i;
  double result = 0.0;

  if (coco_vector_contains_nan(x, number_of_variables))
    return NAN;

  for (i = 0; i < number_of_variables; ++i) {
    double base, exponent, si;

    base = sqrt(alpha);
    exponent = (double)(long)i / ((double)(long)number_of_variables - 1);
    if (best_parameter[i] > 0.0) {
      si = pow(base, exponent);
    } else {
      si = -pow(base, exponent);
    }
    /* boundary handling */
    if (x[i] * best_parameter[i] < 25.0) {
      result += 5.0 * fabs(si) - si * x[i];
    } else {
      result += 5.0 * fabs(si) - si * best_parameter[i];
    }
  }

  return result;
}

/**
 * @brief The linear slope function.
 */
struct f_linear_slope_t final : public coco_problem_s {
  f_linear_slope_t(const size_t number_of_variables, const double *best_parameter)
      : coco_problem_s(number_of_variables, 1, 0) {
    size_t i;
    /* best_parameter will be overwritten below */
    init_from_scalars("linear slope function", -5.0, 5.0, 0.0);
    set_id("%s_d%02lu", "linear_slope", number_of_variables);

    /* Compute best solution */
    for (i = 0; i < number_of_variables; ++i) {
      if (best_parameter[i] < 0.0) {
        this->best_parameter[i] = smallest_values_of_interest[i];
      } else {
        this->best_parameter[i] = largest_values_of_interest[i];
      }
    }
    evaluate_function(this->best_parameter.data(), best_value.data());
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_linear_slope_raw(x, number_of_variables, best_parameter.data());
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the gradient of the linear slope function.
   */
  void evaluate_gradient(const double *x, double *y) override {
    static const double alpha = 100.0;
    double base, exponent, si;
    size_t i;

    (void)x; /* silence (C89) compiliers */
    for (i = 0; i < number_of_variables; ++i) {
      base = sqrt(alpha);
      exponent = (double)(long)i / ((double)(long)number_of_variables - 1);
      if (best_parameter[i] > 0.0) {
        si = pow(base, exponent);
      } else {
        si = -pow(base, exponent);
      }
      y[i] = -si;
    }
  }
};

/**
 * @brief Allocates the basic linear slope problem.
 */
static coco_problem_ptr_t f_linear_slope_allocate(const size_t number_of_variables, const double *best_parameter) {
  return std::make_unique<f_linear_slope_t>(number_of_variables, best_parameter);
}

/**
 * @brief Creates the BBOB linear slope problem.
 */
static coco_problem_ptr_t f_linear_slope_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                               const size_t instance, const long rseed,
                                                               const char *problem_id_template,
                                                               const char *problem_name_template) {
  double fopt;
  coco_problem_ptr_t problem;

  std::vector<double> xopt(dimension);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt.data(), rseed, dimension);
  } else {
    bbob2009_compute_xopt(xopt.data(), rseed, dimension);
  }
  fopt = bbob2009_compute_fopt(function, instance);

  problem = f_linear_slope_allocate(dimension, xopt.data());

  /*if large scale test-bed, normalize by dim*/
  if (coco_strfind(problem_name_template, "BBOB large-scale suite") >= 0) {
    problem = transform_obj_norm_by_dim(std::move(problem));
  }
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("1-separable");

  return problem;
}
