/**
 * @file f_sphere.cpp
 * @brief Implementation of the sphere function and problem.
 */

#include <assert.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_norm_by_dim.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_shift.cpp"

/**
 * @brief Implements the sphere function without connections to any COCO
 * structures.
 */
static double f_sphere_raw(const double *x, const size_t number_of_variables) {

  size_t i = 0;

  double result = 0.0;
  for (i = 0; i < number_of_variables; ++i) {
    result += x[i] * x[i];
  }

  return result;
}

/**
 * @brief The sphere function.
 */
struct f_sphere_t final : public coco_problem_s {
  explicit f_sphere_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("sphere function", -5.0, 5.0, 0.0);
    set_id("%s_d%02zu", "sphere", number_of_variables);

    /* Compute best solution */
    evaluate_function(best_parameter.data(), best_value.data());
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_sphere_raw(x, number_of_variables);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the gradient of the sphere function.
   */
  void evaluate_gradient(const double *x, double *y) override {

    for (size_t i = 0; i < number_of_variables; ++i) {
      y[i] = 2.0 * x[i];
    }
  }
};

/**
 * @brief Allocates the basic sphere problem.
 */
static coco_problem_ptr_t f_sphere_allocate(const size_t number_of_variables) {
  return std::make_unique<f_sphere_t>(number_of_variables);
}

/**
 * @brief Creates the BBOB sphere problem.
 */
static coco_problem_ptr_t f_sphere_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                         const size_t instance, const long rseed,
                                                         const coco_bbob_variant_t variant,
                                                         const char *problem_id_template,
                                                         const char *problem_name_template) {

  coco_problem_ptr_t problem;

  std::vector<double> xopt(dimension);
  bbob_compute_xopt(xopt.data(), rseed, dimension, variant);
  double fopt = bbob2009_compute_fopt(function, instance);

  problem = f_sphere_allocate(dimension);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);

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
