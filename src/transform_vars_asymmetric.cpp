/**
 * @file transform_vars_asymmetric.cpp
 * @brief Implementation of performing an asymmetric transformation on decision values.
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
 * @brief Data type for univariate function tasy_uv.
 */
struct tasy_data {
  double beta;
  size_t i;
  size_t n;
};

/**
 * @brief Univariate asymmetric non-linear transformation.
 */
static double tasy_uv(double xi, tasy_data *d) {
  double yi;
  double exponent;
  if (xi > 0.0) {
    exponent = 1.0 + ((d->beta * (double)(long)d->i) / ((double)(long)d->n - 1.0)) * sqrt(xi);
    yi = pow(xi, exponent);
  } else {
    yi = xi;
  }
  return yi;
}

/**
 * @brief Inverse of asymmetric non-linear transformation tasy_uv obtained with brentq.
 */
static double tasy_uv_inv(double yi, tasy_data *d) {
  double xi;
  xi = brentinv((callback_type)&tasy_uv, yi, d);
  return xi;
}

/**
 * @brief Applies the asymmetric non-linear transformation to the decision values.
 */
struct transform_vars_asymmetric_t final : public coco_problem_transformed_t {
  std::vector<double> asymmetric_x; /**< @brief The transformed decision values */
  const double beta;

  transform_vars_asymmetric_t(coco_problem_ptr_t inner, const double beta)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_asymmetric"),
        asymmetric_x(inner_problem->number_of_variables), beta(beta) {
    if (inner_problem->number_of_objectives > 0 && inner_problem->best_parameter_not_zero()) {
      coco_warning("transform_vars_asymmetric(): 'best_parameter' not updated, set to NAN");
      coco_vector_set_to_nan(inner_problem->best_parameter.data(), inner_problem->number_of_variables);
    }
  }

  /**
   * @brief Multivariate, coordinate-wise, asymmetric non-linear transformation.
   */
  void tasy(const double *x) {
    tasy_data d;

    d.beta = beta;
    d.n = number_of_variables;

    for (size_t i = 0; i < number_of_variables; ++i) {
      d.i = i;
      asymmetric_x[i] = tasy_uv(x[i], &d);
    }
  }

  /**
   * @brief Evaluates the transformed function.
   */
  void evaluate_function(const double *x, double *y) override {
    int feasible;

    if (inner_problem->number_of_objectives == 0) {
      coco_evaluate_function(inner_problem.get(), x, y);
      return;
    }

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(this));
      return;
    }

    tasy(x);

    coco_evaluate_function(inner_problem.get(), asymmetric_x.data(), y);

    if (number_of_constraints > 0) {
      std::vector<double> cons_values(number_of_constraints);
      feasible = is_feasible(x, cons_values.data());
      if (feasible)
        assert(y[0] + 1e-13 >= best_value[0]);
    } else
      assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the transformed constraint.
   */
  void evaluate_constraint(const double *x, double *y, int update_counter) override {
    double exponent;

    if (inner_problem->number_of_constraints == 0) {
      coco_evaluate_constraint_optional_update(inner_problem.get(), x, y, update_counter);
      return;
    }

    if (coco_vector_contains_nan(x, coco_problem_get_dimension(this))) {
      coco_vector_set_to_nan(y, coco_problem_get_number_of_constraints(this));
      return;
    }

    /* FIXME (OME): Old  pre-logger
    tasy(x);

    coco_evaluate_constraint(inner_problem.get(), asymmetric_x, y);
    */
    for (size_t i = 0; i < number_of_variables; ++i) {
      if (x[i] > 0.0) {
        exponent = 1.0 + ((beta * (double)(long)i) / ((double)(long)number_of_variables - 1.0)) * sqrt(x[i]);
        asymmetric_x[i] = pow(x[i], exponent);
      } else {
        asymmetric_x[i] = x[i];
      }
    }
    inner_problem->evaluate_constraint(asymmetric_x.data(), y, update_counter);
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_asymmetric(coco_problem_ptr_t inner_problem, const double beta) {
  return std::make_unique<transform_vars_asymmetric_t>(std::move(inner_problem), beta);
}

/**
 * @brief Returns the asymmetric transformation behind the given problem pointer.
 */
static transform_vars_asymmetric_t *transform_vars_asymmetric_get(coco_problem_t *problem) {
  assert(dynamic_cast<transform_vars_asymmetric_t *>(problem) != nullptr);
  return static_cast<transform_vars_asymmetric_t *>(problem);
}

/**
 * @brief Applies the inverse of the asymmetric transformation tasy to the initial solution.
 *
 *        Takes xopt as input to check the solution remains in the bounds
 *        If not, a curve search is performed
 *        xopt is needed because transform_vars_shift is not yet called
 *        in f_{function}_rotated_c_linear_cons_bbob_problem_allocate
 */
static void transform_inv_initial_asymmetric(coco_problem_t *problem, const double *xopt) {
  size_t j;
  bool is_in_bounds;
  double di;
  double xi;
  double halving_factor = .5;

  tasy_data d;

  std::vector<double> sol(problem->number_of_variables);

  d.beta = transform_vars_asymmetric_get(problem)->beta;
  d.n = problem->number_of_variables;

  j = 0;
  while (1) {

    for (size_t i = 0; i < problem->number_of_variables; ++i) {
      d.i = i;
      di = tasy_uv_inv(problem->initial_solution[i] * pow(halving_factor, (double)(long)j), &d);
      xi = di + xopt[i];
      is_in_bounds = (-5.0 < xi && xi < 5.0);
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
  if (!coco_vector_contains_nan(sol.data(), problem->number_of_variables)) {
    for (size_t i = 0; i < problem->number_of_variables; ++i) {
      problem->initial_solution[i] = sol[i];
    }
  }

}
