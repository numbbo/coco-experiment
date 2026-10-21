/**
 * @file transform_vars_shift.cpp
 * @brief Implementation of shifting all decision values by an offset.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"

/**
 * @brief Shifts all decision values by an offset.
 */
struct transform_vars_shift_t final : public coco_problem_transformed_t {
  const std::vector<double> offset;
  std::vector<double> shifted_x;
  const int shift_constraint_only;

  /**
   * CAVEAT: when shifting the constraint only, the best_value of best_parameter
   *         will get in an inconsistent state.
   */
  transform_vars_shift_t(coco_problem_t *inner_problem, const double *offset, const int shift_constraint_only)
      : coco_problem_transformed_t(inner_problem, "transform_vars_shift"),
        offset(offset, offset + inner_problem->number_of_variables), shifted_x(inner_problem->number_of_variables),
        shift_constraint_only(shift_constraint_only) {
    size_t i;

    /* Update the best parameter */
    for (i = 0; i < number_of_variables; i++)
      best_parameter[i] += this->offset[i];

    /* Update the initial solution if any */
    if (!initial_solution.empty())
      for (i = 0; i < number_of_variables; i++)
        initial_solution[i] += this->offset[i];
  }

  /**
   * @brief Evaluates the transformed objective function.
   */
  void evaluate_function(const double *x, double *y) override {
    size_t i;
    int is_feasible;

    if (inner_problem->number_of_objectives == 0 || shift_constraint_only != 0) {
      coco_evaluate_function(inner_problem, x, y);
      return;
    }

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }

    for (i = 0; i < number_of_variables; ++i) {
      shifted_x[i] = x[i] - offset[i];
    }

    coco_evaluate_function(inner_problem, shifted_x.data(), y);

    if (number_of_constraints > 0) {
      std::vector<double> cons_values(number_of_constraints);
      is_feasible = coco_is_feasible(this, x, cons_values.data());
      if (is_feasible)
        assert(y[0] + 1e-13 >= best_value[0]);
    } else
      assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the transformed constraint function.
   */
  void evaluate_constraint(const double *x, double *y, int update_counter) override {
    size_t i;

    if (inner_problem->number_of_constraints == 0) {
      coco_evaluate_constraint_optional_update(inner_problem, x, y, update_counter);
      return;
    }

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }

    for (i = 0; i < number_of_variables; ++i) {
      shifted_x[i] = x[i] - offset[i];
    }
    inner_problem->evaluate_constraint(shifted_x.data(), y, update_counter);
  }

  /**
   * @brief Evaluates the gradient of the transformed function at x
   */
  void evaluate_gradient(const double *x, double *y) override {
    size_t i;

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }

    for (i = 0; i < number_of_variables; ++i) {
      shifted_x[i] = x[i] - offset[i];
    }
    bbob_evaluate_gradient(inner_problem, shifted_x.data(), y);
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_vars_shift(coco_problem_t *inner_problem, const double *offset,
                                            const int shift_constraint_only) {
  return new transform_vars_shift_t(inner_problem, offset, shift_constraint_only);
}
