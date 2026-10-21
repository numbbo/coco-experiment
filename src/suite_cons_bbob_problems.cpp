/**
 * @file  suite_cons_bbob_problems.cpp
 * @brief Implementation of the problems in the constrained BBOB suite.
 *
 * This suite contains 54 constrained functions in continuous domain
 * which are derived from combining 8 single-objective functions of the
 * noiseless bbob test suite with randomly-generated
 * linear constraints perturbed by nonlinear transformations.
 * Each one of the 8 functions is combined with 6 different numbers of
 * constraints: 1, 2, 6, 6+n/2, 6+n and 6+3n.
 *
 * We consider constrained optimization problems of the form
 *
 *     min f(x) s.t. g_i(x) <= 0, i=1,...,l.
 *
 * The constrained functions are built in a way such that the
 * KKT conditions are satisied at the origin, which is the initial
 * optimal solution. We then translate the constrained function
 * by a random vector in order to move the optimal solution away
 * from the origin.
 *
 * ***************************************************
 * General idea for building the constrained functions
 * ***************************************************
 *
 * 1. Choose a bbob function f whose raw version is pseudoconvex
 *    to be the objective function. (see note below)
 * 2. Remove possible nonlinear transformations from f. (see note below)
 * 3. In order to make sure that the feasible set is not empty,
 *    we choose a direction p that should be feasible.
 *    Define p as the gradient of f at the origin.
 * 4. Define the first constraint g_1(x)=a_1'*x by setting its gradient
 *    to a_1 = -p. By definition, g_1(p) < 0. Thus p is feasible
 *    for g_1.
 * 5. Generate the other constraints randomly while making sure that
 *    p remains feasible for each one.
 * 6. Apply to the whole constrained function the nonlinear transformations
 *    that were removed from the objective function in Step 2. (see note below)
 * 7. Choose a random point xopt and move the optimum away from the
 *    origin to xopt by translating the constrained function by -xopt.
 *
 * @note The pseudoconvexity in Step 1 guarantees that the KKT conditions
 *       are sufficient for optimality.
 *
 * @note The removal of possible nonlinear transformations in Step 2
 *       and posterior application in Step 6 are necessary to make
 *       sure that the KKT conditions are satisfied in the optimum
 *       - until then the origin. As explained in the documentation,
 *       the application of the nonlinear transformations in Step 6
 *       does not affect the location of the optimum.
 *
 * @note a_1 is set to -p, i.e. the negative gradient of f at the origin,
 *       in order to have the KKT conditions easily satisfied.
 *
 * @note Steps 1 and 2 are done within the 'allocate' function of the
 *       objective function, e.g. f_ellipsoid_cons_bbob_problem_allocate().
 *
 * @note Steps 4 and 5 are done within c_linear_cons_bbob_problem_allocate().
 *
 * @note The constrained Rastrigin function's construction differs a bit
 *       from the steps above. Since it is a multimodal function with
 *       well distributed local optima, we choose one of its local optima
 *       to be the global constrained optimum by adding constraints
 *       that pass through that point.
 *
 * @note The testbed provides the user an initial solution which is given
 *       by the feasible direction p scaled by some constant.
 *
 * *************************************************
 * COCO data structure for the constrained functions
 * *************************************************
 *
 * First, we create a coco_problem_t object for the objective function.
 * Then, coco_problem_t objects are created for each constraint and
 * stacked together into a single coco_problem_t object (see "c_linear.cpp").
 * Finally, the coco_problem_t object containing the constraints and
 * the coco_problem_t object containing the objective function are
 * stacked together to form the constrained function.
 *
 */

#include <math.h>

#include "coco.h"
#include "coco_utilities.cpp"
#include "coco_random.cpp"
#include "c_linear.cpp"
#include "f_different_powers.cpp"
#include "f_discus.cpp"
#include "f_ellipsoid.cpp"
#include "f_linear_slope.cpp"
#include "f_rastrigin.cpp"
#include "f_sphere.cpp"

#include "transform_vars_asymmetric.cpp"
#include "transform_vars_oscillate.cpp"

/**
 * @brief Maps the initial solution of the base problem through the inverse of the transformations on top of it.
 *
 * If the result plus xopt is not within (-5, 5) in every coordinate, the initial solution of the base problem
 * is halved until it is (xopt is needed because transform_vars_shift is applied afterwards). The initial
 * solution is left unchanged if the result contains NaN.
 */
static void cons_bbob_inverse_initial_solution(coco_problem_t *problem, const double *xopt) {
  const size_t n = problem->get_number_of_variables();
  const double halving_factor = .5;
  std::vector<double> halved(n), inverted(n), sol(n);
  bool is_in_bounds = true;
  double di = 0.0;

  for (size_t j = 0;; ++j) {
    for (size_t i = 0; i < n; ++i)
      halved[i] = problem->get_initial_solution()[i] * pow(halving_factor, static_cast<double>(j));
    if (!problem->from_base_coordinates(halved.data(), inverted.data()))
      coco_error("cons_bbob_inverse_initial_solution(): the transformations of %s are not invertible",
                 problem->get_id().c_str());
    for (size_t i = 0; i < n; ++i) {
      di = inverted[i];
      const double xi = di + xopt[i];
      is_in_bounds = (-5.0 < xi && xi < 5.0);
      if (!is_in_bounds)
        break;
      sol[i] = di;
    }
    if (is_in_bounds || coco_is_nan(di))
      break;
  }
  if (!coco_vector_contains_nan(sol.data(), n))
    problem->set_initial_solution(sol.data());
}

/**
 * @brief Calculates the objective function type based on the value
 *        of "function"
 */
static size_t obj_function_type(const size_t function) {

  size_t problems_per_obj_function_type = 6;
  return static_cast<size_t>(ceil(static_cast<double>(function) / static_cast<double>(problems_per_obj_function_type)));
}

/**
 * @brief Returns the number of linear constraints associated to the
 *        value of "function"
 */
static size_t nb_of_linear_constraints(const size_t function, const size_t dimension) {

  int problems_per_obj_function_type = 6;

  /* Map "function" value into {1, ..., problems_per_obj_function_type} */
  int p = ((static_cast<int>(function) - 1) % problems_per_obj_function_type) + 1;

  if (p == 1)
    return 1;
  else if (p == 2)
    return 2;
  else if (p == 3)
    return 6;
  else if (p == 4)
    return 6 + dimension / 2;
  else if (p == 5)
    return 6 + dimension;
  else
    return 6 + 3 * dimension;
}

/**
 * @brief Scale feasible direction depending on xopt such that
 *        xopt + feasible_direction remains in [-5, 5].
 *
 *
 */
static void feasible_direction_set_length(double *feasible_direction, const double *xopt, size_t dimension,
                                          long rseed) {
  const long seed_offset = 412;    /* was sampled uniform in 0-999 */
  const double feas_shrink = 0.75; /* scale randomly between 0.75 and 1.0 */
  const double feas_bound = 5.0;

  double r[1], maxabs, maxrel;

  maxabs = maxrel = 0.0;
  for (size_t i = 0; i < dimension; ++i) {
    maxabs = coco_double_max(maxabs, fabs(xopt[i]));
    maxrel = coco_double_max(maxrel, feasible_direction[i] / (feas_bound - xopt[i]));
    maxrel = coco_double_max(maxrel, feasible_direction[i] / (-feas_bound - xopt[i]));
  }
  if (maxabs > 4.01)
    coco_warning("feasible_direction_set_length: a component of fabs(xopt) was greater than 4.01");
  if (maxabs > 5.0)
    coco_error("feasible_direction_set_length: a component of fabs(xopt) was greater than 5.0");
  bbob2009_unif(r, 1, rseed + seed_offset);
  coco_vector_scale(feasible_direction, dimension, feas_shrink + r[0] * (1 - feas_shrink), /* nominator */
                    maxrel);                                                               /* denominator */
}

/**
 * @brief Objective function: sphere
 *        Constraint(s): linear
 */
static coco_problem_ptr_t
f_sphere_c_linear_cons_bbob_problem_allocate(const size_t function, const size_t dimension, const size_t instance,
                                             const size_t number_of_linear_constraints, const long rseed,
                                             double *feasible_direction, const double *xopt,
                                             const char *problem_id_template, const char *problem_name_template) {

  coco_problem_ptr_t problem;
  coco_problem_ptr_t problem_c;

  const std::vector<double> all_zeros(dimension, 0.0);

  /* Create the objective function */
  problem = f_sphere_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                           problem_id_template, problem_name_template);

  problem->evaluate_gradient(all_zeros.data(), feasible_direction);
  feasible_direction_set_length(feasible_direction, xopt, dimension, rseed);

  /* Create the constraints. Use the gradient of the objective
   * function at the origin to build the first constraint.
   */
  problem_c = c_linear_cons_bbob_problem_allocate(function, dimension, instance, number_of_linear_constraints,
                                                  problem_id_template, problem_name_template, feasible_direction);

  const std::string problem_type_temp = problem->get_type();
  const std::string problem_type_constraints = problem_c->get_type();
  problem = coco_problem_stacked_allocate(std::move(problem), std::move(problem_c));

  /* Define problem->best_parameter as the origin and store its
   * objective function value into problem->best_value.
   */
  const std::vector<double> origin(dimension, 0.0);
  problem->set_best_parameter(origin.data());
  problem->evaluate_best_value();

  /* To be sure that everything is correct, reset the f-evaluations
   * and g-evaluations counters to zero.
   */
  problem->reset_evaluation_counters();

  /* Apply a translation to the whole problem so that the constrained
   * minimum is no longer at the origin.
   */
  problem = transform_vars_shift(std::move(problem), xopt, 0);

  /* Construct problem type */
  problem->set_type("%s_%s", problem_type_temp.c_str(), problem_type_constraints.c_str());

  return problem;
}

/**
 * @brief Objective function: ellipsoid
 *        Constraint(s): linear
 */
static coco_problem_ptr_t
f_ellipsoid_c_linear_cons_bbob_problem_allocate(const size_t function, const size_t dimension, const size_t instance,
                                                const size_t number_of_linear_constraints, const long rseed,
                                                double *feasible_direction, const double *xopt,
                                                const char *problem_id_template, const char *problem_name_template) {

  coco_problem_ptr_t problem;
  coco_problem_ptr_t problem_c;

  const std::vector<double> all_zeros(dimension, 0.0);

  /* Create the objective function */
  problem = f_ellipsoid_cons_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                   problem_id_template, problem_name_template);

  problem->evaluate_gradient(all_zeros.data(), feasible_direction);
  feasible_direction_set_length(feasible_direction, xopt, dimension, rseed);

  /* Create the constraints. Use the gradient of the objective
   * function at the origin to build the first constraint.
   */
  problem_c = c_linear_cons_bbob_problem_allocate(function, dimension, instance, number_of_linear_constraints,
                                                  problem_id_template, problem_name_template, feasible_direction);

  const std::string problem_type_temp = problem->get_type();

  /* Build the final constrained function by stacking the objective
   * function and the constraints into one coco_problem_t type.
   */
  const std::string problem_type_constraints = problem_c->get_type();
  problem = coco_problem_stacked_allocate(std::move(problem), std::move(problem_c));

  /* Define problem->best_parameter as the origin and store its
   * objective function value into problem->best_value.
   */
  const std::vector<double> origin(dimension, 0.0);
  problem->set_best_parameter(origin.data());
  problem->evaluate_best_value();

  /* To be sure that everything is correct, reset the f-evaluations
   * and g-evaluations counters to zero.
   */
  problem->reset_evaluation_counters();

  problem = transform_vars_oscillate(std::move(problem));
  cons_bbob_inverse_initial_solution(problem.get(), xopt);

  /* Apply a translation to the whole problem so that the constrained
   * minimum is no longer at the origin.
   */
  problem = transform_vars_shift(std::move(problem), xopt, 0);

  /* Construct problem type */
  problem->set_type("%s_%s", problem_type_temp.c_str(), problem_type_constraints.c_str());

  return problem;
}

/**
 * @brief Objective function: rotated ellipsoid
 *        Constraint(s): linear
 */
static coco_problem_ptr_t f_ellipsoid_rotated_c_linear_cons_bbob_problem_allocate(
    const size_t function, const size_t dimension, const size_t instance, const size_t number_of_linear_constraints,
    const long rseed, double *feasible_direction, const double *xopt, const char *problem_id_template,
    const char *problem_name_template) {

  coco_problem_ptr_t problem;
  coco_problem_ptr_t problem_c;

  const std::vector<double> all_zeros(dimension, 0.0);

  /* Create the objective function */
  problem = f_ellipsoid_rotated_cons_bbob_problem_allocate(
      function, dimension, instance, rseed, coco_bbob_variant_t::bbob, problem_id_template, problem_name_template);

  problem->evaluate_gradient(all_zeros.data(), feasible_direction);
  feasible_direction_set_length(feasible_direction, xopt, dimension, rseed);

  /* Create the constraints. Use the gradient of the objective
   * function at the origin to build the first constraint.
   */
  problem_c = c_linear_cons_bbob_problem_allocate(function, dimension, instance, number_of_linear_constraints,
                                                  problem_id_template, problem_name_template, feasible_direction);

  const std::string problem_type_temp = problem->get_type();

  /* Build the final constrained function by stacking the objective
   * function and the constraints into one coco_problem_t type.
   */
  const std::string problem_type_constraints = problem_c->get_type();
  problem = coco_problem_stacked_allocate(std::move(problem), std::move(problem_c));

  /* Define problem->best_parameter as the origin and store its
   * objective function value into problem->best_value.
   */
  const std::vector<double> origin(dimension, 0.0);
  problem->set_best_parameter(origin.data());
  problem->evaluate_best_value();

  /* To be sure that everything is correct, reset the f-evaluations
   * and g-evaluations counters to zero.
   */
  problem->reset_evaluation_counters();

  problem = transform_vars_oscillate(std::move(problem));
  cons_bbob_inverse_initial_solution(problem.get(), xopt);

  /* Apply a translation to the whole problem so that the constrained
   * minimum is no longer at the origin .
   */
  problem = transform_vars_shift(std::move(problem), xopt, 0);

  /* Construct problem type */
  problem->set_type("%s_%s", problem_type_temp.c_str(), problem_type_constraints.c_str());

  return problem;
}

/**
 * @brief Objective function: linear slope
 *        Constraint(s): linear
 */
static coco_problem_ptr_t
f_linear_slope_c_linear_cons_bbob_problem_allocate(const size_t function, const size_t dimension, const size_t instance,
                                                   const size_t number_of_linear_constraints, const long rseed,
                                                   double *feasible_direction, const double *xopt,
                                                   const char *problem_id_template, const char *problem_name_template) {

  coco_problem_ptr_t problem;
  coco_problem_ptr_t problem_c;

  const std::vector<double> all_zeros(dimension, 0.0);

  /* Create the objective function */
  problem = f_linear_slope_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                 problem_id_template, problem_name_template);

  problem->evaluate_gradient(all_zeros.data(), feasible_direction);
  feasible_direction_set_length(feasible_direction, xopt, dimension, rseed);

  /* Create the constraints. Use the gradient of the objective
   * function at the origin to build the first constraint.
   */
  problem_c = c_linear_cons_bbob_problem_allocate(function, dimension, instance, number_of_linear_constraints,
                                                  problem_id_template, problem_name_template, feasible_direction);

  const std::string problem_type_temp = problem->get_type();

  /* Build the final constrained function by stacking the objective
   * function and the constraints into one coco_problem_t type.
   */
  const std::string problem_type_constraints = problem_c->get_type();
  problem = coco_problem_stacked_allocate(std::move(problem), std::move(problem_c));

  /* Define problem->best_parameter as the origin and store its
   * objective function value into problem->best_value.
   */
  const std::vector<double> origin(dimension, 0.0);
  problem->set_best_parameter(origin.data());

  /* Apply a translation to the whole problem so that the constrained
   * minimum is no longer at the origin.
   */
  problem = transform_vars_shift(std::move(problem), xopt, 1);

  problem->evaluate_best_value();

  /* To be sure that everything is correct, reset the f-evaluations
   * and g-evaluations counters to zero.
   */
  problem->reset_evaluation_counters();

  /* Construct problem type */
  problem->set_type("%s_%s", problem_type_temp.c_str(), problem_type_constraints.c_str());

  return problem;
}

/**
 * @brief Objective function: discus
 *        Constraint(s): linear
 */
static coco_problem_ptr_t
f_discus_c_linear_cons_bbob_problem_allocate(const size_t function, const size_t dimension, const size_t instance,
                                             const size_t number_of_linear_constraints, const long rseed,
                                             double *feasible_direction, const double *xopt,
                                             const char *problem_id_template, const char *problem_name_template) {

  coco_problem_ptr_t problem;
  coco_problem_ptr_t problem_c;

  const std::vector<double> all_zeros(dimension, 0.0);

  /* Create the objective function */
  problem = f_discus_cons_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                problem_id_template, problem_name_template);

  problem->evaluate_gradient(all_zeros.data(), feasible_direction);
  feasible_direction_set_length(feasible_direction, xopt, dimension, rseed);

  /* Create the constraints. Use the gradient of the objective
   * function at the origin to build the first constraint.
   */
  problem_c = c_linear_cons_bbob_problem_allocate(function, dimension, instance, number_of_linear_constraints,
                                                  problem_id_template, problem_name_template, feasible_direction);

  const std::string problem_type_temp = problem->get_type();

  /* Build the final constrained function by stacking the objective
   * function and the constraints into one coco_problem_t type.
   */
  const std::string problem_type_constraints = problem_c->get_type();
  problem = coco_problem_stacked_allocate(std::move(problem), std::move(problem_c));

  /* Define problem->best_parameter as the origin and store its
   * objective function value into problem->best_value
   */
  const std::vector<double> origin(dimension, 0.0);
  problem->set_best_parameter(origin.data());
  problem->evaluate_best_value();

  /* To be sure that everything is correct, reset the f-evaluations
   * and g-evaluations counters to zero.
   */
  problem->reset_evaluation_counters();

  problem = transform_vars_oscillate(std::move(problem));
  cons_bbob_inverse_initial_solution(problem.get(), xopt);

  /* Apply a translation to the whole problem so that the constrained
   * minimum is no longer at the origin.
   */
  problem = transform_vars_shift(std::move(problem), xopt, 0);

  /* Construct problem type */
  problem->set_type("%s_%s", problem_type_temp.c_str(), problem_type_constraints.c_str());

  return problem;
}

/**
 * @brief Objective function: bent cigar
 *        Constraint(s): linear
 */
static coco_problem_ptr_t
f_bent_cigar_c_linear_cons_bbob_problem_allocate(const size_t function, const size_t dimension, const size_t instance,
                                                 const size_t number_of_linear_constraints, const long rseed,
                                                 double *feasible_direction, const double *xopt,
                                                 const char *problem_id_template, const char *problem_name_template) {

  coco_problem_ptr_t problem;
  coco_problem_ptr_t problem_c;

  const std::vector<double> all_zeros(dimension, 0.0);

  /* Create the objective function */
  problem = f_bent_cigar_cons_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                    problem_id_template, problem_name_template);

  problem->evaluate_gradient(all_zeros.data(), feasible_direction);
  feasible_direction_set_length(feasible_direction, xopt, dimension, rseed);

  /* Create the constraints. Use the gradient of the objective
   * function at the origin to build the first constraint.
   */
  problem_c = c_linear_cons_bbob_problem_allocate(function, dimension, instance, number_of_linear_constraints,
                                                  problem_id_template, problem_name_template, feasible_direction);

  const std::string problem_type_temp = problem->get_type();

  /* Build the final constrained function by stacking the objective
   * function and the constraints into one coco_problem_t type.
   */
  const std::string problem_type_constraints = problem_c->get_type();
  problem = coco_problem_stacked_allocate(std::move(problem), std::move(problem_c));

  /* Define problem->best_parameter as the origin and store its
   * objective function value into problem->best_value.
   */
  const std::vector<double> origin(dimension, 0.0);
  problem->set_best_parameter(origin.data());
  problem->evaluate_best_value();

  /* To be sure that everything is correct, reset the f-evaluations
   * and g-evaluations counters to zero.
   */
  problem->reset_evaluation_counters();

  problem = transform_vars_asymmetric(std::move(problem), 0.2);
  cons_bbob_inverse_initial_solution(problem.get(), xopt);

  /* Apply a translation to the whole problem so that the constrained
   * minimum is no longer at the origin.
   */
  problem = transform_vars_shift(std::move(problem), xopt, 0);

  /* Construct problem type */
  problem->set_type("%s_%s", problem_type_temp.c_str(), problem_type_constraints.c_str());

  return problem;
}

/**
 * @brief Objective function: different powers
 *        Constraint(s): linear
 */
static coco_problem_ptr_t f_different_powers_c_linear_cons_bbob_problem_allocate(
    const size_t function, const size_t dimension, const size_t instance, const size_t number_of_linear_constraints,
    const long rseed, double *feasible_direction, const double *xopt, const char *problem_id_template,
    const char *problem_name_template) {

  coco_problem_ptr_t problem;
  coco_problem_ptr_t problem_c;

  const std::vector<double> all_zeros(dimension, 0.0);

  /* Create the objective function */
  problem = f_different_powers_bbob_constrained_problem_allocate(
      function, dimension, instance, rseed, coco_bbob_variant_t::bbob, problem_id_template, problem_name_template);

  problem->evaluate_gradient(all_zeros.data(), feasible_direction);
  feasible_direction_set_length(feasible_direction, xopt, dimension, rseed);

  /* Create the constraints. Use the gradient of the objective
   * function at the origin to build the first constraint.
   */
  problem_c = c_linear_cons_bbob_problem_allocate(function, dimension, instance, number_of_linear_constraints,
                                                  problem_id_template, problem_name_template, feasible_direction);

  const std::string problem_type_temp = problem->get_type();

  /* Build the final constrained function by stacking the objective
   * function and the constraints into one coco_problem_t type.
   */
  const std::string problem_type_constraints = problem_c->get_type();
  problem = coco_problem_stacked_allocate(std::move(problem), std::move(problem_c));

  /* Define problem->best_parameter as the origin and store its
   * objective function value into problem->best_value.
   */
  const std::vector<double> origin(dimension, 0.0);
  problem->set_best_parameter(origin.data());
  problem->evaluate_best_value();

  /* To be sure that everything is correct, reset the f-evaluations
   * and g-evaluations counters to zero.
   */
  problem->reset_evaluation_counters();

  /* Apply a translation to the whole problem so that the constrained
   * minimum is no longer at the origin.
   */
  problem = transform_vars_shift(std::move(problem), xopt, 0);

  /* Construct problem type */
  problem->set_type("%s_%s", problem_type_temp.c_str(), problem_type_constraints.c_str());

  return problem;
}

/**
 * @brief Objective function: rastrigin
 *        Constraint(s): linear
 */
static coco_problem_ptr_t
f_rastrigin_c_linear_cons_bbob_problem_allocate(const size_t function, const size_t dimension, const size_t instance,
                                                const size_t number_of_linear_constraints, const long rseed,
                                                double *feasible_direction, const double *xopt,
                                                const char *problem_id_template, const char *problem_name_template) {

  coco_problem_ptr_t problem;
  coco_problem_ptr_t problem_c;

  const std::vector<double> all_zeros(dimension, 0.0);

  /* Create the objective function */
  problem = f_rastrigin_cons_bbob_problem_allocate(function, dimension, instance, rseed, problem_id_template,
                                                   problem_name_template);

  problem->evaluate_gradient(all_zeros.data(), feasible_direction);
  feasible_direction_set_length(feasible_direction, xopt, dimension, rseed);

  /* Create the constraints. Use the feasible direction above
   * to build the first constraint.
   */
  problem_c = c_linear_cons_bbob_problem_allocate(function, dimension, instance, number_of_linear_constraints,
                                                  problem_id_template, problem_name_template, feasible_direction);

  const std::string problem_type_temp = problem->get_type();

  /* Build the final constrained function by stacking the objective
   * function and the constraints into one coco_problem_t type.
   */
  const std::string problem_type_constraints = problem_c->get_type();
  problem = coco_problem_stacked_allocate(std::move(problem), std::move(problem_c));

  /* Define problem->best_parameter as the origin and store its
   * objective function value into problem->best_value.
   */
  const std::vector<double> origin(dimension, 0.0);
  problem->set_best_parameter(origin.data());
  problem->evaluate_best_value();

  /* To be sure that everything is correct, reset the f-evaluations
   * and g-evaluations counters to zero.
   */
  problem->reset_evaluation_counters();

  problem = transform_vars_asymmetric(std::move(problem), 0.2);
  problem = transform_vars_oscillate(std::move(problem));
  cons_bbob_inverse_initial_solution(problem.get(), xopt);

  /* Apply a translation to the whole problem so that the constrained
   * minimum is no longer at the origin.
   */
  problem = transform_vars_shift(std::move(problem), xopt, 0);

  /* Construct problem type */
  problem->set_type("%s_%s", problem_type_temp.c_str(), problem_type_constraints.c_str());

  return problem;
}

/**
 * @brief Objective function: rastrigin
 *        Constraint(s): linear
 */
static coco_problem_ptr_t f_rastrigin_rotated_c_linear_cons_bbob_problem_allocate(
    const size_t function, const size_t dimension, const size_t instance, const size_t number_of_linear_constraints,
    const long rseed, double *feasible_direction, const double *xopt, const char *problem_id_template,
    const char *problem_name_template) {

  coco_problem_ptr_t problem;
  coco_problem_ptr_t problem_c;

  const std::vector<double> all_zeros(dimension, 0.0);

  /* Create the objective function */
  problem = f_rastrigin_rotated_cons_bbob_problem_allocate(function, dimension, instance, rseed, problem_id_template,
                                                           problem_name_template);

  problem->evaluate_gradient(all_zeros.data(), feasible_direction);
  feasible_direction_set_length(feasible_direction, xopt, dimension, rseed);

  /* Create the constraints. Use the feasible direction above
   * to build the first constraint.
   */
  problem_c = c_linear_cons_bbob_problem_allocate(function, dimension, instance, number_of_linear_constraints,
                                                  problem_id_template, problem_name_template, feasible_direction);

  const std::string problem_type_temp = problem->get_type();

  /* Build the final constrained function by stacking the objective
   * function and the constraints into one coco_problem_t type.
   */
  const std::string problem_type_constraints = problem_c->get_type();
  problem = coco_problem_stacked_allocate(std::move(problem), std::move(problem_c));

  /* Define problem->best_parameter as the origin and store its
   * objective function value into problem->best_value.
   */
  const std::vector<double> origin(dimension, 0.0);
  problem->set_best_parameter(origin.data());
  problem->evaluate_best_value();

  /* To be sure that everything is correct, reset the f-evaluations
   * and g-evaluations counters to zero.
   */
  problem->reset_evaluation_counters();

  problem = transform_vars_asymmetric(std::move(problem), 0.2);
  problem = transform_vars_oscillate(std::move(problem));
  cons_bbob_inverse_initial_solution(problem.get(), xopt);

  /* Apply a translation to the whole problem so that the constrained
   * minimum is no longer at the origin.
   */
  problem = transform_vars_shift(std::move(problem), xopt, 0);

  /* Construct problem type */
  problem->set_type("%s_%s", problem_type_temp.c_str(), problem_type_constraints.c_str());

  return problem;
}
