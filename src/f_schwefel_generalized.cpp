/**
 * @file f_schwefel_generalized.cpp
 * @brief Implementation of the Schwefel function and problem.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_norm_by_dim.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_affine.cpp"
#include "transform_vars_conditioning.cpp"
#include "transform_vars_scale.cpp"
#include "transform_vars_shift.cpp"
#include "transform_vars_x_hat.cpp"
#include "transform_vars_z_hat.cpp"

/**
 * @brief Implements the Schwefel function for large scale.
 */
static double f_schwefel_generalized_raw(const double *x,
                                         const size_t number_of_variables) {

  size_t i = 0;
  double penalty, sum;
  const double Schwefel_constant = 4.189828872724339;

  /* Boundary handling*/
  penalty = 0.0;
  for (i = 0; i < number_of_variables; ++i) {
    const double tmp = fabs(x[i]) - 500.0;
    if (tmp > 0.0)
      penalty += tmp * tmp;
  }

  /* Computation core */
  sum = 0.0;
  for (i = 0; i < number_of_variables; ++i) {
    sum += x[i] * sin(sqrt(fabs(x[i])));
  }
  double result = 0.01 * (penalty + Schwefel_constant * 100. - sum / static_cast<double>(number_of_variables));

  return result;
}

/**
 * @brief The Schwefel function.
 */
struct f_schwefel_generalized_t final : public coco_problem_s {
  explicit f_schwefel_generalized_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("Schwefel function", -5.0, 5.0, 420.96874633);
    set_id("%s_d%02zu", "schwefel", number_of_variables);

    /* Compute best solution: best_parameter[i] = 200 * fabs(xopt[i]) */
    evaluate_function(best_parameter.data(), best_value.data());
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_schwefel_generalized_raw(x, number_of_variables);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Allocates the basic Schwefel problem.
 */
static coco_problem_ptr_t f_schwefel_generalized_allocate(const size_t number_of_variables) {
  return std::make_unique<f_schwefel_generalized_t>(number_of_variables);
}

/**
 * @brief Creates the BBOB Schwefel problem.
 */
static coco_problem_ptr_t f_schwefel_generalized_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                       const size_t instance, const long rseed,
                                                                       const char *problem_id_template,
                                                                       const char *problem_name_template) {
  coco_problem_ptr_t problem;
  std::vector<double> tmp1(dimension);
  std::vector<double> tmp2(dimension);

  std::vector<double> xopt(dimension);
  double fopt = bbob2009_compute_fopt(function, instance);
  bbob2009_unif(tmp1.data(), dimension, rseed);
  for (size_t i = 0; i < dimension; ++i) {
    xopt[i] = 0.5 * 4.2096874637;
    if (tmp1[i] - 0.5 < 0) {
      xopt[i] *= -1;
    }
  }

  for (size_t i = 0; i < dimension; ++i) {
    tmp1[i] = -2 * fabs(xopt[i]);
    tmp2[i] = 2 * fabs(xopt[i]);
  }

  problem = f_schwefel_generalized_allocate(dimension);
  problem = transform_vars_scale(std::move(problem), 100);
  problem = transform_vars_shift(std::move(problem), tmp1.data(), 0);
  problem = transform_vars_conditioning(std::move(problem), 10.0);
  problem = transform_vars_shift(std::move(problem), tmp2.data(), 0);
  problem = transform_vars_z_hat(std::move(problem), xopt.data());
  problem = transform_vars_scale(std::move(problem), 2);
  problem = transform_vars_x_hat(std::move(problem), rseed);
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("5-weakly-structured");

  return problem;
}
