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

/**
 * @brief Applies the affine transformation x |-> Mx + b to the decision values.
 */
class transform_vars_affine_t final : public coco_problem_transformed_t {
public:
  transform_vars_affine_t(coco_problem_ptr_t inner, coco_matrix_t M, std::vector<double> b)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_affine"), M(std::move(M)), b(std::move(b)),
        transformed_x(inner_problem->get_number_of_variables()) {
    /* Update the best parameter by computing
       problem->best_parameter = M^T * (inner_problem->best_parameter - b).

       The update takes place only if the best parameter or b are different than zero
       and the transformation matrix M is orthogonal.
    */
    if (inner_problem->best_parameter_not_zero() ||
        !coco_vector_is_zero(this->b.data(), inner_problem->get_number_of_variables())) {
      if (!this->M.is_orthogonal())
        coco_warning("transform_vars_affine(): rotation matrix is not orthogonal. Best parameter not updated");
      else {
        for (size_t i = 0; i < inner_problem->get_number_of_variables(); ++i) {
          transformed_x[i] = inner_problem->get_best_parameter()[i] - this->b[i];
        }
        for (size_t i = 0; i < this->number_of_variables; ++i) {
          best_parameter[i] = 0;
          for (size_t j = 0; j < inner_problem->get_number_of_variables(); ++j) {
            best_parameter[i] += this->M(j, i) * transformed_x[j];
          }
        }
      }
    }
  }

  /**
   * @brief Computes Mx + b and stores it in transformed_x.
   */
  void transform(const double *x) { M.multiply_add(x, b.data(), transformed_x.data()); }

  /**
   * @brief Evaluates the transformed objective function.
   */
  void evaluate_function(const double *x, double *y) override {
    if (inner_problem->get_number_of_objectives() == 0) {
      inner_problem->evaluate_function(x, y);
      return;
    }

    transform(x);
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

    transform(x);
    inner_problem->evaluate_constraint(transformed_x.data(), y, update_counter);
  }

  /**
   * @brief Evaluates the gradient of the transformed function.
   */
  void evaluate_gradient(const double *x, double *y) override {

    std::vector<double> gradient(inner_problem->get_number_of_variables(), 0.0);

    transform(x);
    inner_problem->evaluate_gradient(transformed_x.data(), y);

    /* grad_(f o g )(x), where g(x) = M * x + b, equals to
     * M^T * grad_f(M *x + b)
     */
    for (size_t j = 0; j < inner_problem->get_number_of_variables(); ++j) {
      for (size_t i = 0; i < inner_problem->get_number_of_variables(); ++i) {
        gradient[j] += y[i] * M(i, j);
      }
    }
    for (size_t i = 0; i < inner_problem->get_number_of_variables(); ++i)
      y[i] = gradient[i];
  }

  /** @brief Computes M^T (inner_x - b), which is the inverse only if M is orthogonal. */
  bool inverse(const double *inner_x, double *x) const override {
    if (!M.is_orthogonal())
      return false;
    std::vector<double> shifted_x(inner_problem->get_number_of_variables());
    for (size_t i = 0; i < inner_problem->get_number_of_variables(); ++i)
      shifted_x[i] = inner_x[i] - b[i];
    M.multiply_transposed(shifted_x.data(), x);
    return true;
  }

private:
  const coco_matrix_t M; /**< @brief Has number_of_variables columns and inner number of variables rows */
  const std::vector<double> b;
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
