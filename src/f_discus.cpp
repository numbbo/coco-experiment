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
  size_t i;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
    return NAN;

  result = condition * x[0] * x[0];
  for (i = 1; i < number_of_variables; ++i) {
    result += x[i] * x[i];
  }

  return result;
}

/**
 * @brief The discus function.
 */
struct f_discus_t final : public coco_problem_s {
  explicit f_discus_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    coco_problem_init_from_scalars(this, "discus function", -5.0, 5.0, 0.0);
    coco_problem_set_id(this, "%s_d%02lu", "discus", number_of_variables);

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
    size_t i;

    y[0] = condition * 2.0 * x[0];
    for (i = 1; i < number_of_variables; ++i)
      y[i] = 2.0 * x[i];
  }
};

/**
 * @brief Allocates the basic discus problem.
 */
static coco_problem_t *f_discus_allocate(const size_t number_of_variables) {
  return new f_discus_t(number_of_variables);
}

/**
 * @brief Creates the BBOB discus problem.
 */
static coco_problem_t *f_discus_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                      const size_t instance, const long rseed,
                                                      const char *problem_id_template,
                                                      const char *problem_name_template) {

  double *xopt, fopt;
  coco_problem_t *problem = NULL;

  double *M = coco_allocate_vector(dimension * dimension);
  double *b = coco_allocate_vector(dimension);
  double **rot1;

  xopt = coco_allocate_vector(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt, rseed, dimension);
  } else {
    bbob2009_compute_xopt(xopt, rseed, dimension);
  }

  rot1 = bbob2009_allocate_matrix(dimension, dimension);
  bbob2009_compute_rotation(rot1, rseed + 1000000, dimension);
  bbob2009_copy_rotation_matrix(rot1, M, b, dimension);
  bbob2009_free_matrix(rot1, dimension);

  problem = f_discus_allocate(dimension);
  problem = transform_vars_oscillate(problem);
  problem = transform_vars_affine(problem, M, b, dimension);
  problem = transform_vars_shift(problem, xopt, 0);
  problem = transform_obj_shift(problem, fopt);

  coco_problem_set_id(problem, problem_id_template, function, instance,
                      dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance,
                        dimension);
  coco_problem_set_type(problem, "3-ill-conditioned");

  coco_free_memory(M);
  coco_free_memory(b);
  coco_free_memory(xopt);
  return problem;
}

/**
 * @brief Creates the discus problem for the constrained BBOB suite.
 */
static coco_problem_t *f_discus_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                           const size_t instance, const long rseed,
                                                           const char *problem_id_template,
                                                           const char *problem_name_template) {

  double *xopt, fopt;
  coco_problem_t *problem = NULL;

  double *M = coco_allocate_vector(dimension * dimension);
  double *b = coco_allocate_vector(dimension);
  double **rot1;

  xopt = coco_allocate_vector(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt, rseed, dimension);
  } else {
    bbob2009_compute_xopt(xopt, rseed, dimension);
  }

  rot1 = bbob2009_allocate_matrix(dimension, dimension);
  bbob2009_compute_rotation(rot1, rseed + 1000000, dimension);
  bbob2009_copy_rotation_matrix(rot1, M, b, dimension);
  bbob2009_free_matrix(rot1, dimension);

  problem = f_discus_allocate(dimension);
  problem = transform_vars_affine(problem, M, b, dimension);
  problem = transform_vars_shift(problem, xopt, 0);
  problem = transform_obj_shift(problem, fopt);

  coco_problem_set_id(problem, problem_id_template, function, instance,
                      dimension);
  coco_problem_set_name(problem, problem_name_template, function, instance,
                        dimension);
  coco_problem_set_type(problem, "3-ill-conditioned");

  coco_free_memory(M);
  coco_free_memory(b);
  coco_free_memory(xopt);
  return problem;
}
