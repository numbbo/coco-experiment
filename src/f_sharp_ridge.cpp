/**
 * @file f_sharp_ridge.cpp
 * @brief Implementation of the sharp ridge function and problem.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_affine.cpp"
#include "transform_vars_shift.cpp"

/**
 * @brief Implements the sharp ridge function without connections to any COCO
 * structures.
 */
static double f_sharp_ridge_raw(const double *x,
                                const size_t number_of_variables) {

  static const double alpha = 100.0;
  const double d_vars_40 = 1.0; /* generalized: number_of_variables <= 40 ? 1 :
                                   number_of_variables / 40.0; */
  const size_t vars_40 = coco_double_to_size_t(ceil(d_vars_40));
  size_t i = 0;

  assert(number_of_variables > 1);

  double result = 0.0;
  for (i = vars_40; i < number_of_variables; ++i) {
    result += x[i] * x[i];
  }
  result = alpha * sqrt(result / d_vars_40);
  for (i = 0; i < vars_40; ++i)
    result += x[i] * x[i] / d_vars_40;

  return result;
}

/**
 * @brief The sharp ridge function.
 */
struct f_sharp_ridge_t final : public coco_problem_s {
  explicit f_sharp_ridge_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("sharp ridge function", -5.0, 5.0, 0.0);
    set_id("%s_d%02zu", "sharp_ridge", number_of_variables);

    /* Compute best solution */
    evaluate_function(best_parameter.data(), best_value.data());
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_sharp_ridge_raw(x, number_of_variables);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Allocates the basic sharp ridge problem.
 */
static coco_problem_ptr_t f_sharp_ridge_allocate(const size_t number_of_variables) {
  return std::make_unique<f_sharp_ridge_t>(number_of_variables);
}

/**
 * @brief Creates the BBOB sharp ridge problem.
 */
static coco_problem_ptr_t f_sharp_ridge_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                              const size_t instance, const long rseed,
                                                              const coco_bbob_variant_t variant,
                                                              const char *problem_id_template,
                                                              const char *problem_name_template) {
  coco_problem_ptr_t problem;
  coco_matrix_t M(dimension, dimension);
  std::vector<double> b(dimension);

  std::vector<double> xopt(dimension);
  double fopt = bbob2009_compute_fopt(function, instance);
  bbob_compute_xopt(xopt.data(), rseed, dimension, variant);

  const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
  const coco_matrix_t rot2 = bbob2009_compute_rotation(rseed, dimension);
  std::vector<double> scales(dimension);
  for (size_t k = 0; k < dimension; ++k) {
    const double exponent = 1.0 * static_cast<int>(k) / (static_cast<double>(dimension) - 1.0);
    scales[k] = pow(sqrt(10), exponent);
  }
  M = coco_matrix_product_with_diagonal(rot1, scales, rot2);
  problem = f_sharp_ridge_allocate(dimension);
  problem = transform_obj_shift(std::move(problem), fopt);
  problem = transform_vars_affine(std::move(problem), M, b);
  problem = transform_vars_shift(std::move(problem), xopt.data(), 0);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("3-ill-conditioned");

  return problem;
}
