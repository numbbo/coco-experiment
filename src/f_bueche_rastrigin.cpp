/**
 * @file f_bueche_rastrigin.cpp
 * @brief Implementation of the Bueche-Rastrigin function and problem.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_norm_by_dim.cpp"
#include "transform_obj_penalize.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_brs.cpp"
#include "transform_vars_oscillate.cpp"
#include "transform_vars_shift.cpp"

/**
 * @brief Implements the Bueche-Rastrigin function without connections to any
 * COCO structures.
 */
static double f_bueche_rastrigin_raw(const double *x,
                                     const size_t number_of_variables) {
  double tmp = 0., tmp2 = 0.;
  size_t i;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables)) {
    return NAN;
  }

  for (i = 0; i < number_of_variables; ++i) {
    tmp += cos(2 * coco_pi * x[i]);
    tmp2 += x[i] * x[i];
  }
  result = 10.0 * ((double)(long)number_of_variables - tmp) + tmp2 + 0;
  return result;
}

/**
 * @brief The Bueche-Rastrigin function.
 */
struct f_bueche_rastrigin_t final : public coco_problem_s {
  explicit f_bueche_rastrigin_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    coco_problem_init_from_scalars(this, "Bueche-Rastrigin function", -5.0, 5.0, 0.0);
    coco_problem_set_id(this, "%s_d%02lu", "bueche-rastrigin", number_of_variables);

    /* Compute best solution */
    evaluate_function(best_parameter.data(), best_value.data());
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_bueche_rastrigin_raw(x, number_of_variables);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Allocates the basic Bueche-Rastrigin problem.
 */
static coco_problem_ptr_t f_bueche_rastrigin_allocate(const size_t number_of_variables) {
  return std::make_unique<f_bueche_rastrigin_t>(number_of_variables);
}

/**
 * @brief Creates the BBOB Bueche-Rastrigin problem.
 */
static coco_problem_ptr_t f_bueche_rastrigin_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                   const size_t instance, const long rseed,
                                                                   const char *problem_id_template,
                                                                   const char *problem_name_template) {
  double fopt;
  coco_problem_ptr_t problem;

  const double penalty_factor = 100.0;
  size_t i;

  std::vector<double> xopt(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt.data(), rseed, dimension);
  } else {
    bbob2009_compute_xopt(xopt.data(), rseed, dimension);
  }

  /* OME: This step is in the legacy C code but _not_ in the function
   * description. */
  for (i = 0; i < dimension; i += 2) {
    xopt[i] = fabs(xopt[i]);
  }

  problem = f_bueche_rastrigin_allocate(dimension);
  problem = transform_vars_brs(std::move(problem));
  problem = transform_vars_oscillate(std::move(problem));
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);

  /*if large scale test-bed, normalize by dim*/
  if (coco_strfind(problem_name_template, "BBOB large-scale suite") >= 0) {
    problem = transform_obj_norm_by_dim(std::move(problem));
  }
  problem = transform_obj_shift(std::move(problem), fopt);
  problem = transform_obj_penalize(std::move(problem), penalty_factor);

  coco_problem_set_id(problem.get(), problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem.get(), problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem.get(), "1-separable");

  return problem;
}
