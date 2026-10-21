/**
 * @file transform_vars_oscillate.cpp
 * @brief Implementation of oscillating the decision values.
 * @author ??
 * @author Paul Dufosse
 * @note Edited to fulfill needs from the constrained test bed.
 */

#include <math.h>
#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "brentq.cpp"

/**
 * @brief Data type for univariate function tosz_uv
 */
typedef struct {
  double alpha;
} tosz_data;

/**
 * @brief Univariate oscillating non-linear transformation.
 */
static double tosz_uv(double xi, tosz_data *d) {
  double yi;
  double tmp, base;
  if (xi > 0.0) {
    tmp = log(xi) / d->alpha;
    base = exp(tmp + 0.49 * (sin(tmp) + sin(0.79 * tmp)));
    yi = pow(base, d->alpha);
  } else if (xi < 0.0) {
    tmp = log(-xi) / d->alpha;
    base = exp(tmp + 0.49 * (sin(0.55 * tmp) + sin(0.31 * tmp)));
    yi = -pow(base, d->alpha);
  } else {
    yi = 0.0;
  }
  return yi;
}

/**
 * @brief Inverse of oscillating non-linear transformation tosz_uv obtained with brentq.
 */
static double tosz_uv_inv(double yi, tosz_data *d) {
  double xi;
  xi = brentinv((callback_type)&tosz_uv, yi, d);
  return xi;
}

/**
 * @brief Applies the oscillating non-linear transformation to the decision values.
 */
struct transform_vars_oscillate_t final : public coco_problem_transformed_t {
  const double alpha;
  std::vector<double> oscillated_x;

  explicit transform_vars_oscillate_t(coco_problem_t *inner_problem)
      : coco_problem_transformed_t(inner_problem, "transform_vars_oscillate"), alpha(0.1),
        oscillated_x(inner_problem->number_of_variables) {}

  /**
   * @brief Multivariate, coordinate-wise, oscillating non-linear transformation.
   */
  void tosz(const double *x) {
    size_t i;
    tosz_data d;

    d.alpha = alpha;

    for (i = 0; i < number_of_variables; ++i) {
      oscillated_x[i] = tosz_uv(x[i], &d);
    }
  }

  /**
   * @brief Evaluates the transformed objective functions.
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

    tosz(x);

    coco_evaluate_function(inner_problem, oscillated_x.data(), y);

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
   * @brief Evaluates the transformed constraints.
   */
  void evaluate_constraint(const double *x, double *y, int update_counter) override {
    double tmp, base;
    size_t i;

    if (inner_problem->number_of_constraints == 0) {
      coco_evaluate_constraint_optional_update(inner_problem, x, y, update_counter);
      return;
    }

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_constraints(this));
      return;
    }

    for (i = 0; i < number_of_variables; ++i) {
      if (x[i] > 0.0) {
        tmp = log(x[i]) / alpha;
        base = exp(tmp + 0.49 * (sin(tmp) + sin(0.79 * tmp)));
        oscillated_x[i] = pow(base, alpha);
      } else if (x[i] < 0.0) {
        tmp = log(-x[i]) / alpha;
        base = exp(tmp + 0.49 * (sin(0.55 * tmp) + sin(0.31 * tmp)));
        oscillated_x[i] = -pow(base, alpha);
      } else {
        oscillated_x[i] = 0.0;
      }
    }
    inner_problem->evaluate_constraint(oscillated_x.data(), y, update_counter);
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_vars_oscillate(coco_problem_t *inner_problem) {
  return new transform_vars_oscillate_t(inner_problem);
}

/**
 * @brief Returns the oscillating transformation behind the given problem pointer.
 */
static transform_vars_oscillate_t *transform_vars_oscillate_get(coco_problem_t *problem) {
  assert(dynamic_cast<transform_vars_oscillate_t *>(problem) != NULL);
  return static_cast<transform_vars_oscillate_t *>(problem);
}

/**
 * @brief Applies the inverse of the oscillating transformation tasy(tosz(.)) to the initial solution.
 *
 *        Takes xopt as input to check the solution remains in the bounds
 *        If not, a curve search is performed
 *        xopt is needed because transform_vars_shift is not yet called
 *        in f_{function}_rotated_c_linear_cons_bbob_problem_allocate
 */
static void transform_inv_initial_oscillate(coco_problem_t *problem, const double *xopt) {
  size_t i;
  size_t j;
  int is_in_bounds;
  double di = 0.0;
  double xi;
  double *sol = NULL;
  double halving_factor = .5;

  tosz_data *d;

  sol = coco_allocate_vector(problem->number_of_variables);
  d = (tosz_data *)coco_allocate_memory(sizeof(*d));

  d->alpha = transform_vars_oscillate_get(problem)->alpha;

  j = 0;
  while (1) {

    for (i = 0; i < problem->number_of_variables; ++i) {
      di = tosz_uv_inv(problem->initial_solution[i] * pow(halving_factor, (double)(long)j), d);
      xi = di + xopt[i];
      is_in_bounds = (int)(-5.0 < xi && xi < 5.0);
      /* Line search for the inverse-transformed feasible initial solution
         to remain within the bounds
        */
      if (!is_in_bounds) {
        j = j + 1;
        break;
      }
      sol[i] = di;
    }
    if (!is_in_bounds && !coco_is_nan(di)) {
      continue;
    } else {
      break;
    }
  }
  if (!coco_vector_contains_nan(sol, problem->number_of_variables)) {
    for (i = 0; i < problem->number_of_variables; ++i) {
      problem->initial_solution[i] = sol[i];
    }
  }
  coco_free_memory(d);
  coco_free_memory(sol);
}
