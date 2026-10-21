/**
 * @file transform_vars_discretize.cpp
 *
 * @brief Implementation of transforming a continuous problem to a mixed-integer problem by making some
 * of its variables discrete. The integer variables are considered as bounded (any variable outside the
 * decision space is mapped to the closest boundary point), while the continuous ones are treated as
 * unbounded.
 *
 * @note The first problem->number_of_integer_variables are integer, while the rest are continuous.
 *
 * The discretization works as follows. Consider the case where the interval [l, u] of the inner problem
 * needs to be discretized to n integer values of the outer problem. First, [l, u] is discretized to n
 * integers by placing the integers so that there is a (u-l)/(n+1) distance between them (and the border
 * points). Then, the transformation is shifted so that the optimum aligns with the closest integer. In
 * this way, we make sure that if the optimum is within [l, u], so are all the shifted points.
 *
 * When evaluating such a problem, the x values of the integer variables are first discretized. Any value
 * x < 0 is mapped to 0 and any value x > (n-1) is mapped to (n-1).
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"

/**
 * @brief Makes the first number_of_integer_variables variables of the inner problem discrete.
 */
struct transform_vars_discretize_t final : public coco_problem_transformed_t {
  std::vector<double> offset;

  transform_vars_discretize_t(coco_problem_t *inner_problem, const double *smallest_values_of_interest,
                              const double *largest_values_of_interest, const size_t number_of_integer_variables)
      : coco_problem_transformed_t(inner_problem, "transform_vars_discretize"),
        offset(inner_problem->number_of_variables) {
    double l, u, inner_l, inner_u, outer_l, outer_u;
    double outer_xopt, inner_xopt, inner_approx_xopt;
    const double precision_offset = 1e-7; /* Needed to avoid issues with rounding doubles */
    int n;
    size_t i;

    assert(number_of_integer_variables > 0);
    this->number_of_integer_variables = number_of_integer_variables;

    for (i = 0; i < number_of_variables; i++) {
      assert(smallest_values_of_interest[i] < largest_values_of_interest[i]);
      this->smallest_values_of_interest[i] = smallest_values_of_interest[i];
      this->largest_values_of_interest[i] = largest_values_of_interest[i];
      offset[i] = 0;
      if (i < number_of_integer_variables) {
        /* Compute the offset for integer variables */
        outer_l = this->smallest_values_of_interest[i];
        outer_u = this->largest_values_of_interest[i];
        l = inner_problem->smallest_values_of_interest[i];
        u = inner_problem->largest_values_of_interest[i];
        n = coco_double_to_int(outer_u) - coco_double_to_int(outer_l) + 1; /* number of integer values */
        assert(n > 1);
        inner_l = l + (u - l) / (n + 1);
        inner_u = u - (u - l) / (n + 1);
        /* Find the location of the optimum in the coordinates of the outer problem */
        inner_xopt = inner_problem->best_parameter[i];
        outer_xopt = outer_l + (outer_u - outer_l) * (inner_xopt - inner_l) / (inner_u - inner_l);
        outer_xopt = coco_double_round(outer_xopt + precision_offset);
        /* Make sure you the bounds are respected */
        if (outer_xopt < outer_l)
          outer_xopt = outer_l;
        if (outer_xopt > outer_u)
          outer_xopt = outer_u;
        best_parameter[i] = outer_xopt;
        /* Find the value corresponding to outer_xopt in the coordinates of the inner problem */
        inner_approx_xopt = inner_l + (inner_u - inner_l) * (outer_xopt - outer_l) / (outer_u - outer_l);
        /* Compute the difference between the inner_approx_xopt and inner_xopt */
        offset[i] = inner_approx_xopt - inner_xopt;
      }
    }

    if (number_of_constraints > 0)
      coco_error("transform_vars_discretize(): Constraints not supported yet.");
  }

  /**
   * @brief Evaluates the transformed objective function.
   */
  void evaluate_function(const double *x, double *y) override {
    size_t i;
    double l, u, inner_l, inner_u, outer_l, outer_u;
    int n;

    if (inner_problem->number_of_objectives == 0) {
      coco_evaluate_function(inner_problem, x, y);
      return;
    }

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }

    /* Transform x to fit in the discretized space */
    std::vector<double> discretized_x(x, x + number_of_variables);
    for (i = 0; i < number_of_integer_variables; ++i) {
      outer_l = smallest_values_of_interest[i];
      outer_u = largest_values_of_interest[i];
      l = inner_problem->smallest_values_of_interest[i];
      u = inner_problem->largest_values_of_interest[i];
      n = coco_double_to_int(outer_u) - coco_double_to_int(outer_l) + 1; /* number of integer values in this coordinate */
      assert(n > 1);
      inner_l = l + (u - l) / (n + 1);
      inner_u = u - (u - l) / (n + 1);
      /* Make sure you the bounds are respected */
      discretized_x[i] = coco_double_round(x[i]);
      if (discretized_x[i] < outer_l)
        discretized_x[i] = outer_l;
      if (discretized_x[i] > outer_u)
        discretized_x[i] = outer_u;
      discretized_x[i] = inner_l + (inner_u - inner_l) * (discretized_x[i] - outer_l) / (outer_u - outer_l) - offset[i];
    }

    coco_evaluate_function(inner_problem, discretized_x.data(), y);
  }

  /** @brief Constraints are not supported yet (TODO?). */
  void evaluate_constraint(const double *x, double *y, int update_counter) override {
    coco_problem_s::evaluate_constraint(x, y, update_counter);
  }

  /** @brief Gradients are not supported yet (TODO?). */
  void evaluate_gradient(const double *x, double *y) override { coco_problem_s::evaluate_gradient(x, y); }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_vars_discretize(coco_problem_t *inner_problem,
                                                 const double *smallest_values_of_interest,
                                                 const double *largest_values_of_interest,
                                                 const size_t number_of_integer_variables) {
  return new transform_vars_discretize_t(inner_problem, smallest_values_of_interest, largest_values_of_interest,
                                         number_of_integer_variables);
}
