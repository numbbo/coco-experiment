/**
 * @file f_discus.cpp
 * @brief Implementation of the discus function and problem.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_affine.cpp"
#include "transform_vars_oscillate.cpp"
#include "transform_vars_shift.cpp"

/**
 * @brief Implements the discus function without connections to any COCO
 * structures.
 */
static double f_discus_raw(const double *x, const size_t number_of_variables) {

  static const double condition = 1.0e6;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
    return NAN;

  result = condition * x[0] * x[0];
  for (size_t i = 1; i < number_of_variables; ++i) {
    result += x[i] * x[i];
  }

  return result;
}

/**
 * @brief The discus function.
 */
struct f_discus_t final : public coco_problem_s {
  explicit f_discus_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("discus function", -5.0, 5.0, 0.0);
    set_id("%s_d%02lu", "discus", number_of_variables);

    /* Compute best solution */
    evaluate_function(best_parameter.data(), best_value.data());
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_discus_raw(x, number_of_variables);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the gradient of the discus function.
   */
  void evaluate_gradient(const double *x, double *y) override {

    static const double condition = 1.0e6;

    y[0] = condition * 2.0 * x[0];
    for (size_t i = 1; i < number_of_variables; ++i)
      y[i] = 2.0 * x[i];
  }
};

/**
 * @brief Allocates the basic discus problem.
 */
static coco_problem_ptr_t f_discus_allocate(const size_t number_of_variables) {
  return std::make_unique<f_discus_t>(number_of_variables);
}

/**
 * @brief Creates the BBOB discus problem.
 */
static coco_problem_ptr_t f_discus_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                         const size_t instance, const long rseed,
                                                         const coco_bbob_variant_t variant,
                                                         const char *problem_id_template,
                                                         const char *problem_name_template) {

  double fopt;
  coco_problem_ptr_t problem;

  coco_matrix_t M(dimension, dimension);
  std::vector<double> b(dimension);

  std::vector<double> xopt(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  bbob_compute_xopt(xopt.data(), rseed, dimension, variant);

  const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
  M = rot1;

  problem = f_discus_allocate(dimension);
  problem = transform_vars_oscillate(std::move(problem));
  problem = transform_vars_affine(std::move(problem), M, b);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("3-ill-conditioned");

  return problem;
}

/**
 * @brief Creates the discus problem for the constrained BBOB suite.
 */
static coco_problem_ptr_t f_discus_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                              const size_t instance, const long rseed,
                                                              const coco_bbob_variant_t variant,
                                                              const char *problem_id_template,
                                                              const char *problem_name_template) {

  double fopt;
  coco_problem_ptr_t problem;

  coco_matrix_t M(dimension, dimension);
  std::vector<double> b(dimension);

  std::vector<double> xopt(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  bbob_compute_xopt(xopt.data(), rseed, dimension, variant);

  const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
  M = rot1;

  problem = f_discus_allocate(dimension);
  problem = transform_vars_affine(std::move(problem), M, b);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("3-ill-conditioned");

  return problem;
}
