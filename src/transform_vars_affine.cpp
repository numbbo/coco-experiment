/**
 * @file transform_vars_affine.cpp
 * @brief Implementation of performing an affine transformation on decision values.
 *
 * x |-> Mx + b
 *
 * Currently, the best parameter is transformed correctly only in the simple
 * cases where M is orthogonal which is always the case for the `bbob`
 * functions. How to code this for general transformations of the above form,
 * see https://github.com/numbbo/coco/issues/814#issuecomment-303724400
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "transform_functions.cpp"

/**
 * @brief Applies the affine transformation x |-> Mx + b to the decision values.
 */
class transform_vars_affine_t final : public coco_problem_transformed_t {
public:
  transform_vars_affine_t(coco_problem_ptr_t inner, coco_matrix_t M, std::vector<double> b)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_affine"), affine(std::move(M), std::move(b)),
        transformed_x(inner_problem->get_number_of_variables()) {
    update_best_parameter();
  }

  /**
   * @brief Evaluates the transformed objective function.
   */
  void evaluate_function(const double *x, double *y) override {
    if (inner_problem->get_number_of_objectives() == 0) {
      inner_problem->evaluate_function(x, y);
      return;
    }

    affine.apply(x, transformed_x.data());
    inner_problem->evaluate_function(transformed_x.data(), y);

    assert_not_below_best_value(x, y);
  }

  /**
   * @brief Evaluates the transformed constraint.
   */
  void evaluate_constraint(const double *x, double *y, int update_counter) override {
    if (inner_problem->get_number_of_constraints() == 0) {
      inner_problem->evaluate_constraint(x, y, update_counter);
      return;
    }

    affine.apply(x, transformed_x.data());
    inner_problem->evaluate_constraint(transformed_x.data(), y, update_counter);
  }

  /**
   * @brief Evaluates the gradient of the transformed function.
   */
  void evaluate_gradient(const double *x, double *y) override {

    std::vector<double> gradient(inner_problem->get_number_of_variables(), 0.0);

    affine.apply(x, transformed_x.data());
    inner_problem->evaluate_gradient(transformed_x.data(), y);

    /* grad_(f o g )(x), where g(x) = M * x + b, equals to
     * M^T * grad_f(M *x + b)
     */
    const coco_matrix_t &M = affine.matrix();
    for (size_t j = 0; j < inner_problem->get_number_of_variables(); ++j) {
      for (size_t i = 0; i < inner_problem->get_number_of_variables(); ++i) {
        gradient[j] += y[i] * M(i, j);
      }
    }
    for (size_t i = 0; i < inner_problem->get_number_of_variables(); ++i)
      y[i] = gradient[i];
  }

  /**
   * @brief Solves M x = inner_x - b, which is M^T (inner_x - b) if M is orthogonal. Not possible if M is singular.
   */
  bool inverse(const double *inner_x, double *x) const override { return affine.inverse(inner_x, x); }

private:
  const vars_affine_t affine;
  std::vector<double> transformed_x; /**< @brief The transformed decision values */
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_affine(coco_problem_ptr_t inner_problem, coco_matrix_t M,
                                                std::vector<double> b) {
  /*
   * TODOs:
   * - Calculate new smallest/largest values of interest?
   * - Resize bounds vectors if input and output dimensions do not match
   */
  assert(M.number_of_rows() == inner_problem->get_number_of_variables() && b.size() == M.number_of_rows());
  return std::make_unique<transform_vars_affine_t>(std::move(inner_problem), std::move(M), std::move(b));
}
