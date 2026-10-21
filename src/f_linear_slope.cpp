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
static double f_linear_slope_raw(const double *x, const size_t number_of_variables, const double *best_parameter,
                                 const double *slopes) {
  double result = 0.0;

  for (size_t i = 0; i < number_of_variables; ++i) {
    const double si = best_parameter[i] > 0.0 ? slopes[i] : -slopes[i];
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
class f_linear_slope_t final : public coco_problem_s {
public:
  f_linear_slope_t(const size_t number_of_variables, const double *best_parameter)
      : coco_problem_s(number_of_variables, 1, 0), slopes(number_of_variables) {
    static const double alpha = 100.0;
    for (size_t i = 0; i < number_of_variables; ++i)
      slopes[i] = pow(sqrt(alpha), static_cast<double>(i) / (static_cast<double>(number_of_variables) - 1));

    /* best_parameter will be overwritten below */
    init_from_scalars("linear slope function", -5.0, 5.0, 0.0);
    set_id("%s_d%02zu", "linear_slope", number_of_variables);

    /* Compute best solution */
    for (size_t i = 0; i < number_of_variables; ++i) {
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
    y[0] = f_linear_slope_raw(x, number_of_variables, best_parameter.data(), slopes.data());
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the gradient of the linear slope function.
   */
  void evaluate_gradient(const double *, double *y) override {
    for (size_t i = 0; i < number_of_variables; ++i)
      y[i] = best_parameter[i] > 0.0 ? -slopes[i] : slopes[i];
  }

private:
  std::vector<double> slopes; /**< @brief |s_i| = sqrt(alpha)^(i / (n - 1)) */
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
                                                               const coco_bbob_variant_t variant,
                                                               const char *problem_id_template,
                                                               const char *problem_name_template) {
  coco_problem_ptr_t problem;

  std::vector<double> xopt(dimension);
  bbob_compute_xopt(xopt.data(), rseed, dimension, variant);
  double fopt = bbob2009_compute_fopt(function, instance);

  problem = f_linear_slope_allocate(dimension, xopt.data());

  /*if large scale test-bed, normalize by dim*/
  if (variant == coco_bbob_variant_t::largescale) {
    problem = transform_obj_norm_by_dim(std::move(problem));
  }
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("1-separable");

  return problem;
}
