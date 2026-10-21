/**
 * @file transform_vars_affine.cpp
 * @brief Implementation of performing an affine transformation on decision values.
 *
 * x |-> Mx + b <br>
 * The matrix M is stored in row-major format.
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
struct transform_vars_affine_t final : public coco_problem_transformed_t {
  const std::vector<double> M; /**< @brief Has number_of_variables columns and inner number of variables rows */
  const std::vector<double> b;
  std::vector<double> transformed_x; /**< @brief The transformed decision values */

  transform_vars_affine_t(coco_problem_t *inner_problem, const double *M, const double *b,
                          const size_t number_of_variables)
      : coco_problem_transformed_t(inner_problem, "transform_vars_affine"),
        M(M, M + inner_problem->number_of_variables * number_of_variables),
        b(b, b + inner_problem->number_of_variables), transformed_x(inner_problem->number_of_variables) {
    size_t i, j;
    /* Update the best parameter by computing
       problem->best_parameter = M^T * (inner_problem->best_parameter - b).

       The update takes place only if the best parameter or b are different than zero
       and the transformation matrix M is orthogonal.
    */
    if (coco_problem_best_parameter_not_zero(inner_problem) ||
        !coco_vector_is_zero(this->b.data(), inner_problem->number_of_variables)) {
      if (!coco_is_orthogonal(this->M.data(), this->number_of_variables, inner_problem->number_of_variables))
        coco_warning("transform_vars_affine(): rotation matrix is not orthogonal. Best parameter not updated");
      else {
        for (i = 0; i < inner_problem->number_of_variables; ++i) {
          transformed_x[i] = inner_problem->best_parameter[i] - this->b[i];
        }
        for (i = 0; i < this->number_of_variables; ++i) {
          best_parameter[i] = 0;
          for (j = 0; j < inner_problem->number_of_variables; ++j) {
            best_parameter[i] += this->M[j * this->number_of_variables + i] * transformed_x[j];
          }
        }
      }
    }
  }

  /**
   * @brief Computes Mx + b and stores it in transformed_x.
   */
  void transform(const double *x) {
    size_t i, j;
    for (i = 0; i < inner_problem->number_of_variables; ++i) {
      const double *current_row = M.data() + i * number_of_variables;
      transformed_x[i] = b[i];
      for (j = 0; j < number_of_variables; ++j) {
        transformed_x[i] += x[j] * current_row[j];
      }
    }
  }

  /**
   * @brief Evaluates the transformed objective function.
   */
  void evaluate_function(const double *x, double *y) override {
    double *cons_values;
    int is_feasible;

    if (inner_problem->number_of_objectives == 0) {
      coco_evaluate_function(inner_problem, x, y);
      return;
    }

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }

    transform(x);
    coco_evaluate_function(inner_problem, transformed_x.data(), y);

    if (number_of_constraints > 0) {
      cons_values = coco_allocate_vector(number_of_constraints);
      is_feasible = coco_is_feasible(this, x, cons_values);
      coco_free_memory(cons_values);
      if (is_feasible)
        assert(y[0] + 1e-13 >= best_value[0]);
    } else
      assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the transformed constraint.
   */
  void evaluate_constraint(const double *x, double *y, int update_counter) override {
    if (inner_problem->number_of_constraints == 0) {
      coco_evaluate_constraint_optional_update(inner_problem, x, y, update_counter);
      return;
    }

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_constraints(this));
      return;
    }

    transform(x);
    inner_problem->evaluate_constraint(transformed_x.data(), y, update_counter);
  }

  /**
   * @brief Evaluates the gradient of the transformed function.
   */
  void evaluate_gradient(const double *x, double *y) override {
    size_t i, j;
    const double *current_row;
    double *gradient;

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }

    gradient = coco_allocate_vector(inner_problem->number_of_variables);
    for (i = 0; i < inner_problem->number_of_variables; ++i)
      gradient[i] = 0.0;

    transform(x);
    bbob_evaluate_gradient(inner_problem, transformed_x.data(), y);

    /* grad_(f o g )(x), where g(x) = M * x + b, equals to
     * M^T * grad_f(M *x + b)
     */
    for (j = 0; j < inner_problem->number_of_variables; ++j) {
      for (i = 0; i < inner_problem->number_of_variables; ++i) {
        current_row = M.data() + i * number_of_variables;
        gradient[j] += y[i] * current_row[j];
      }
    }
    for (i = 0; i < inner_problem->number_of_variables; ++i)
      y[i] = gradient[i];
    coco_free_memory(gradient);
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_vars_affine(coco_problem_t *inner_problem, const double *M, const double *b,
                                             const size_t number_of_variables) {
  /*
   * TODOs:
   * - Calculate new smallest/largest values of interest?
   * - Resize bounds vectors if input and output dimensions do not match
   */
  return new transform_vars_affine_t(inner_problem, M, b, number_of_variables);
}
