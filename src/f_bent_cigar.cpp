/**
 * @file f_bent_cigar.cpp
 * @brief Implementation of the bent cigar function and problem.
 */

#include <assert.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_affine.cpp"
#include "transform_vars_asymmetric.cpp"
#include "transform_vars_shift.cpp"

/**
 * @brief Implements the bent cigar function without connections to any COCO
 * structures.
 */
static double f_bent_cigar_raw(const double *x,
                               const size_t number_of_variables) {

  static const double condition = 1.0e6;
  size_t i;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
    return NAN;

  result = x[0] * x[0];
  for (i = 1; i < number_of_variables; ++i) {
    result += condition * x[i] * x[i];
  }
  return result;
}

/**
 * @brief The bent cigar function.
 */
struct f_bent_cigar_t final : public coco_problem_s {
  explicit f_bent_cigar_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("bent cigar function", -5.0, 5.0, 0.0);
    set_id("%s_d%02lu", "bent_cigar", number_of_variables);

    /* Compute best solution */
    evaluate_function(best_parameter.data(), best_value.data());
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_bent_cigar_raw(x, number_of_variables);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the gradient of the bent cigar function.
   */
  void evaluate_gradient(const double *x, double *y) override {

    static const double condition = 1.0e6;
    size_t i;

    y[0] = 2.0 * x[0];
    for (i = 1; i < number_of_variables; ++i)
      y[i] = 2.0 * condition * x[i];
  }
};

/**
 * @brief Allocates the basic bent cigar problem.
 */
static coco_problem_ptr_t f_bent_cigar_allocate(const size_t number_of_variables) {
  return std::make_unique<f_bent_cigar_t>(number_of_variables);
}

/**
 * @brief Creates the BBOB bent cigar problem.
 */
static coco_problem_ptr_t f_bent_cigar_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                             const size_t instance, const long rseed,
                                                             const char *problem_id_template,
                                                             const char *problem_name_template) {

  double fopt;
  coco_problem_ptr_t problem;

  std::vector<double> M(dimension * dimension);
  std::vector<double> b(dimension);

  std::vector<double> xopt(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt.data(), rseed + 1000000, dimension);
  } else {
    bbob2009_compute_xopt(xopt.data(), rseed + 1000000, dimension);
  }

  coco_matrix_t rot1(dimension, dimension);
  bbob2009_compute_rotation(rot1.rows(), rseed + 1000000, dimension);
  bbob2009_copy_rotation_matrix(rot1.rows(), M.data(), b.data(), dimension);

  problem = f_bent_cigar_allocate(dimension);
  problem = transform_obj_shift(std::move(problem), fopt);
  problem = transform_vars_affine(std::move(problem), M.data(), b.data(), dimension);
  problem = transform_vars_asymmetric(std::move(problem), 0.5);
  problem = transform_vars_affine(std::move(problem), M.data(), b.data(), dimension);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("3-ill-conditioned");

  return problem;
}

/**
 * @brief Creates the bent cigar problem for the constrained BBOB suite.
 */
static coco_problem_ptr_t f_bent_cigar_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                  const size_t instance, const long rseed,
                                                                  const char *problem_id_template,
                                                                  const char *problem_name_template) {

  double fopt;
  coco_problem_ptr_t problem;

  std::vector<double> M(dimension * dimension);
  std::vector<double> b(dimension);

  std::vector<double> xopt(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(xopt.data(), rseed + 1000000, dimension);
  } else {
    bbob2009_compute_xopt(xopt.data(), rseed + 1000000, dimension);
  }

  coco_matrix_t rot1(dimension, dimension);
  bbob2009_compute_rotation(rot1.rows(), rseed + 1000000, dimension);
  bbob2009_copy_rotation_matrix(rot1.rows(), M.data(), b.data(), dimension);

  problem = f_bent_cigar_allocate(dimension);
  problem = transform_obj_shift(std::move(problem), fopt);
  problem = transform_vars_affine(std::move(problem), M.data(), b.data(), dimension);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("3-ill-conditioned");

  return problem;
}
