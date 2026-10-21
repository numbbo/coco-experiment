/**
 * @file coco_problem.cpp
 * @brief Definitions of functions regarding COCO problems.
 */

#include <float.h>
#include "coco.h"
#include "coco_internal.h"

#include "coco_utilities.cpp"

/***********************************************************************************************************/

/**
 * @name Methods regarding the basic COCO problem
 */
/**@{*/
/**
 * Evaluates the problem function, increases the number of evaluations and updates the best observed value
 * and the best observed evaluation number.
 *
 * @note Both x and y must point to correctly sized allocated memory regions.
 *
 * @param problem The given COCO problem.
 * @param x The decision vector.
 * @param y The objective vector that is the result of the evaluation (in single-objective problems only the
 * first vector item is being set).
 */
void coco_evaluate_function(coco_problem_t *problem, const double *x, double *y) {
  /* implements a safer version of problem->evaluate(problem, x, y) */
  assert(problem != nullptr);
  const std::span<const double> x_span(x, coco_problem_get_dimension(problem));
  const std::span<double> y_span(y, coco_problem_get_number_of_objectives(problem));

  /* Set objective vector to INFINITY if the decision vector contains any INFINITY values */
  for (const double xi : x_span) {
    if (coco_is_inf(xi)) {
      for (double &yj : y_span) {
        yj = fabs(xi);
      }
      return;
    }
  }

  /* Set objective vector to NAN if the decision vector contains any NAN values */
  if (coco_vector_contains_nan(x_span)) {
    coco_vector_set_to_nan(y_span);
    return;
  }

  problem->evaluate_function(x_span, y_span);
  problem->count_evaluation(x_span, y_span);
}

/**
 * Evaluates the problem constraint.
 *
 * @note Both x and y must point to correctly sized allocated memory regions.
 *
 * @param problem The given COCO problem.
 * @param x The decision vector.
 * @param y The vector of constraints that is the result of the evaluation.
 */
void coco_evaluate_constraint_optional_update(coco_problem_t *problem, const double *x, double *y, int update_counter) {
  /* implements a safer version of problem->evaluate(problem, x, y) */
  assert(problem != nullptr);
  const std::span<const double> x_span(x, coco_problem_get_dimension(problem));
  const std::span<double> y_span(y, coco_problem_get_number_of_constraints(problem));

  /* Set constraints vector to INFINITY if the decision vector contains any INFINITY values */
  for (const double xi : x_span) {
    if (coco_is_inf(xi)) {
      for (double &yj : y_span) {
        yj = fabs(xi);
      }
      return;
    }
  }

  /* Set constraints vector to NAN if the decision vector contains any NAN values */
  if (coco_vector_contains_nan(x_span)) {
    coco_vector_set_to_nan(y_span);
    return;
  }

  problem->evaluate_constraint(x_span, y_span, update_counter);
  if (update_counter)
    problem->count_constraint_evaluation();
}

void coco_evaluate_constraint(coco_problem_t *problem, const double *x, double *y) {
  coco_evaluate_constraint_optional_update(problem, x, y, 1);
}

/**
 * Evaluates and logs the given solution (as the coco_evaluate_function), but does not return the evaluated
 * value.
 *
 * @note x must point to a correctly sized allocated memory region.

 * @param problem The given COCO problem.
 * @param x The decision vector.
 */
void coco_recommend_solution(coco_problem_t *problem, const double *x) {
  assert(problem != nullptr);
  problem->recommend_solution({x, coco_problem_get_dimension(problem)});
}

/**
 * @name Default implementations of the problem methods
 *
 * They report that the method is not implemented for the problem.
 */
/**@{*/
void coco_problem_s::evaluate_function(std::span<const double>, std::span<double>) {
  coco_error("coco_evaluate_function(): No function implemented for problem %s", problem_id.c_str());
}

void coco_problem_s::evaluate_constraint(std::span<const double>, std::span<double>, int) {
  coco_error("coco_evaluate_constraint_optional_update(): No constraint function implemented for problem %s",
             problem_id.c_str());
}

void coco_problem_s::evaluate_gradient(std::span<const double>, std::span<double>) {
  coco_error("evaluate_gradient(): No gradient function implemented for problem %s", problem_id.c_str());
}

void coco_problem_s::recommend_solution(std::span<const double>) {
  coco_warning("coco_recommend_solution(): No recommend solution function implemented for problem %s",
               problem_id.c_str());
}

void coco_problem_s::assert_not_below_best_value([[maybe_unused]] std::span<const double> x,
                                                 [[maybe_unused]] std::span<const double> y) {
#ifndef NDEBUG
  if (number_of_constraints == 0 || is_feasible(x, {}))
    assert(y[0] + 1e-13 >= best_value[0]);
#endif
}

bool coco_problem_s::from_base_coordinates(std::span<const double> base_x, std::span<double> x) const {
  assert(base_x.size() == x.size());
  std::copy(base_x.begin(), base_x.end(), x.begin());
  return true;
}
/**@}*/
/***********************************************************************************************************/

/**
 * @brief Constructs a problem for the given number of variables, number of objectives and number of constraints.
 */
coco_problem_s::coco_problem_s(const size_t number_of_variables, const size_t number_of_objectives,
                               const size_t number_of_constraints) {
  /* Initialize fields to sane/safe defaults */
  this->number_of_variables = number_of_variables;
  this->number_of_objectives = number_of_objectives;
  this->number_of_constraints = number_of_constraints;
  smallest_values_of_interest.assign(number_of_variables, 0.0);
  largest_values_of_interest.assign(number_of_variables, 0.0);
  number_of_integer_variables = 0; /* No integer variables by default */

  noisy = false;
  if (number_of_objectives > 0)
    last_noise_free_values.assign(number_of_objectives, 0.0);
  if (number_of_objectives > 1) {
    opt_known = false; /* Optimum of multi-objective problems is unknown by default */
    best_value.assign(number_of_objectives, 0.0);
    nadir_value.assign(number_of_objectives, 0.0);
  } else {
    opt_known = true; /* Optimum of single-objective problems is known by default */
    best_parameter.assign(number_of_variables, 0.0);
    best_value.assign(1, 0.0);
  }
  tainted = false;

  evaluations = 0;
  evaluations_constraints = 0;
  final_target_delta = 1e-8; /* in case to be modified by the benchmark */
  best_observed_fvalue = DBL_MAX;
  best_observed_evaluation = 0;
  suite = nullptr; /* To be initialized in the coco_suite_get_problem_from_indices() function */
  suite_dep_index = 0;
  suite_dep_function = 0;
  suite_dep_instance = 0;
}

/**
 * @brief Copies the properties of the 'other' problem, i.e. everything except its behaviour, data and noise.
 *
 * Both problems must have the same number of variables, objectives and constraints.
 */
void coco_problem_s::copy_properties_from(const coco_problem_s &other) {
  smallest_values_of_interest = other.smallest_values_of_interest;
  largest_values_of_interest = other.largest_values_of_interest;
  if (!other.best_parameter.empty())
    best_parameter = other.best_parameter;
  tainted = other.tainted;

  number_of_integer_variables = other.number_of_integer_variables;

  initial_solution = other.initial_solution;

  opt_known = other.opt_known;
  if (!other.best_value.empty())
    best_value = other.best_value;
  if (!other.nadir_value.empty())
    nadir_value = other.nadir_value;

  problem_name = other.problem_name;
  problem_id = other.problem_id;
  problem_type = other.problem_type;

  evaluations = other.evaluations;
  evaluations_constraints = other.evaluations_constraints;
  final_target_delta = other.final_target_delta;
  best_observed_fvalue = other.best_observed_fvalue;
  best_observed_evaluation = other.best_observed_evaluation;

  suite = other.suite;
  suite_dep_index = other.suite_dep_index;
  suite_dep_function = other.suite_dep_function;
  suite_dep_instance = other.suite_dep_instance;
}

void coco_problem_s::count_evaluation(std::span<const double> x, std::span<const double> y) {
  evaluations++; /* each derived class has its own counter, only the most outer will be visible */

  /* A little bit of bookkeeping */
  if (y[0] < best_observed_fvalue) {
    bool feasible = true;
    if (number_of_constraints > 0) {
      std::vector<double> z(number_of_constraints);
      feasible = is_feasible(x, z);
    }
    if (feasible) {
      best_observed_fvalue = y[0];
      best_observed_evaluation = evaluations;
    }
  }
}

void coco_problem_s::set_suite_dependent(coco_suite_t *suite, const size_t index, const size_t function,
                                         const size_t instance) {
  this->suite = suite;
  suite_dep_index = index;
  suite_dep_function = function;
  suite_dep_instance = instance;
}

/**
 * @brief Initializes a single-objective problem using scalar values for smallest_value_of_interest,
 * largest_value_of_interest and best_parameter. Assumes all variables are continuous.
 */
void coco_problem_s::init_from_scalars(std::string_view name, const double smallest, const double largest,
                                       const double best) {

  assert(number_of_objectives == 1);
  assert(number_of_constraints == 0);
  problem_name = name;

  for (size_t i = 0; i < number_of_variables; ++i) {
    smallest_values_of_interest[i] = smallest;
    largest_values_of_interest[i] = largest;
    best_parameter[i] = best;
  }
  number_of_integer_variables = 0;
}

void coco_problem_free(coco_problem_t *problem) {
  assert(problem != nullptr);
  delete problem;
}

/***********************************************************************************************************/

/**
 * @brief Checks whether the given string is in the right format to be a problem_id.
 *
 * No non-alphanumeric characters besides '-' and '_' are allowed.
 */
static bool coco_problem_id_is_fine(const std::string &id) {
  for (const char c : id) {
    if (('A' <= c) && (c <= 'Z'))
      continue;
    if (('a' <= c) && (c <= 'z'))
      continue;
    if ((c == '_') || (c == '-'))
      continue;
    if (('0' <= c) && (c <= '9'))
      continue;
    return false;
  }
  return true;
}

/**
 * @brief Sets the problem_id using formatted printing (as in printf).
 *
 * Verifies that the problem_id is in the correct format.
 */
void coco_problem_s::set_id(const char *format, ...) {
  va_list args;

  va_start(args, format);
  problem_id = coco_string_vformat(format, args);
  va_end(args);
  if (!coco_problem_id_is_fine(problem_id)) {
    coco_error("Problem id should only contain standard chars, not like '%s'", problem_id.c_str());
  }
}

/**
 * @brief Sets the problem_name using formatted printing (as in printf).
 */
void coco_problem_s::set_name(const char *format, ...) {
  va_list args;

  va_start(args, format);
  problem_name = coco_string_vformat(format, args);
  va_end(args);
}

/**
 * @brief Sets the problem_type using formatted printing (as in printf).
 */
void coco_problem_s::set_type(const char *format, ...) {
  va_list args;

  va_start(args, format);
  problem_type = coco_string_vformat(format, args);
  va_end(args);
}

size_t coco_problem_get_evaluations(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return problem->get_evaluations();
}

size_t coco_problem_get_evaluations_constraints(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return problem->get_evaluations_constraints();
}

/**
 * @note Can be used to prevent unnecessary burning of CPU time.
 */
int coco_problem_final_target_hit(const coco_problem_t *problem) {
  assert(problem != nullptr);
  if (coco_problem_get_number_of_objectives(problem) != 1 || coco_problem_get_evaluations(problem) < 1)
    return 0;
  if (problem->get_best_value().empty())
    return 0;
  return problem->get_best_observed_fvalue() <= problem->get_best_value()[0] + problem->get_final_target_delta() ? 1
                                                                                                                 : 0;
}
/**
 * @note Tentative...
 */
double coco_problem_get_best_observed_fvalue1(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return problem->get_best_observed_fvalue();
}

/**
 * @brief Returns the optimal function value of the problem
 */
double coco_problem_get_best_value(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(!problem->get_best_value().empty());
  return problem->get_best_value()[0];
}

/**
 * @note This function breaks the black-box property: the returned  value is not
 * meant to be used by the optimization algorithm.
 */
double coco_problem_get_final_target_fvalue1(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(!problem->get_best_value().empty());
  return problem->get_best_value()[0] + problem->get_final_target_delta();
}

/**
 * @note Do not modify the returned string! If you free the problem, the returned pointer becomes invalid.
 * When in doubt, copy the returned value.
 */
const char *coco_problem_get_name(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(!problem->get_name().empty());
  return problem->get_name().c_str();
}

/**
 * The ID is guaranteed to contain only characters in the set [a-z0-9_-]. It should therefore be safe to use
 * it to construct filenames or other identifiers.
 *
 * Each problem ID should be unique within each benchmark suite.
 *
 * @note Do not modify the returned string! If you free the problem, the returned pointer becomes invalid.
 * When in doubt, copy the returned value.
 */
const char *coco_problem_get_id(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(!problem->get_id().empty());
  return problem->get_id().c_str();
}

const char *coco_problem_get_type(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(!problem->get_type().empty());
  return problem->get_type().c_str();
}

size_t coco_problem_get_dimension(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(problem->get_number_of_variables() > 0);
  return problem->get_number_of_variables();
}

size_t coco_problem_get_number_of_objectives(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return problem->get_number_of_objectives();
}

size_t coco_problem_get_number_of_constraints(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return problem->get_number_of_constraints();
}

const double *coco_problem_get_smallest_values_of_interest(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(!problem->get_smallest_values_of_interest().empty());
  return problem->get_smallest_values_of_interest().data();
}

const double *coco_problem_get_largest_values_of_interest(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(!problem->get_largest_values_of_interest().empty());
  return problem->get_largest_values_of_interest().data();
}

size_t coco_problem_get_number_of_integer_variables(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return problem->get_number_of_integer_variables();
}

const double *coco_problem_get_largest_fvalues_of_interest(const coco_problem_t *problem) {
  assert(problem != nullptr);
  if (problem->get_number_of_objectives() == 1)
    coco_error(
        "coco_problem_get_largest_fvalues_of_interest(): f-values of interest undefined for single-objective problems");
  if (problem->get_nadir_value().empty())
    coco_error("coco_problem_get_largest_fvalues_of_interest(): f-values of interest undefined");
  return problem->get_nadir_value().data();
}

/**
 * Copies problem->initial_solution into initial_solution if not null,
 * otherwise the center of the problem's region of interest is the
 * initial solution. Takes care of rounding the solution in case of integer variables.
 *
 * @param problem The given COCO problem.
 * @param initial_solution The pointer to the initial solution being set by this method.
 */
void coco_problem_get_initial_solution(const coco_problem_t *problem, double *initial_solution) {


  assert(problem != nullptr);
  if (!problem->get_initial_solution().empty()) {
    for (size_t i = 0; i < problem->get_number_of_variables(); ++i)
      initial_solution[i] = problem->get_initial_solution()[i];
  } else {
    assert(!problem->get_smallest_values_of_interest().empty());
    assert(!problem->get_largest_values_of_interest().empty());
    for (size_t i = 0; i < problem->get_number_of_variables(); ++i)
      initial_solution[i] =
          problem->get_smallest_values_of_interest()[i] +
          0.5 * (problem->get_largest_values_of_interest()[i] - problem->get_smallest_values_of_interest()[i]);
    if (problem->get_number_of_integer_variables() > 0) {
      for (size_t i = 0; i < problem->get_number_of_integer_variables(); ++i) {
        initial_solution[i] = coco_double_round(initial_solution[i]);
      }
    }
  }
}

size_t coco_problem_get_suite_dep_index(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return problem->get_suite_dep_index();
}

static size_t coco_problem_get_suite_dep_function(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(problem->get_suite_dep_function() > 0);
  return problem->get_suite_dep_function();
}

static size_t coco_problem_get_suite_dep_instance(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(problem->get_suite_dep_instance() > 0);
  return problem->get_suite_dep_instance();
}
/**@}*/

/**
 * Copies problem->best_parameter into best_parameter and returns 1 if not null,
 * otherwise returns 0. Only works for single objective problems.
 *
 * @param problem The given COCO problem.
 * @param best_parameter The pointer to a vector into which the best parameter is copied.
 *
 * @returns true if successful, otherwise false.
 */
int coco_problem_get_best_parameter(coco_problem_t *problem, double *best_parameter) {
  assert(problem != nullptr);

  if (problem->get_best_parameter().empty() || problem->get_number_of_objectives() != 1) {
    return 0;
  }

  // No more benchmarking on this problem!
  problem->taint();

  for (size_t i = 0; i < problem->get_number_of_variables(); ++i) {
    best_parameter[i] = problem->get_best_parameter()[i];
  }

  return 1;
}

/**
 * Returns whether a problem is tainted, i. e. if it's best parameter has been retrieved.
 *
 * @param problem The given COCO problem.
 *
 * @returns true if problem is tainted, otherwise returns false.re
 */
int coco_problem_is_tainted(const coco_problem_t *problem) { return problem->is_tainted(); }

void bbob_problem_best_parameter_print(const coco_problem_t *problem) {
  coco_file_t file;
  assert(problem != nullptr);
  assert(!problem->get_best_parameter().empty());
  if (file.open("._bbob_problem_best_parameter.txt", "w")) {
    for (size_t i = 0; i < problem->get_number_of_variables(); ++i)
      fprintf(file, " %.16f ", problem->get_best_parameter()[i]);
  }
}

/* We need a forward declaration of this for `bbob_biobj_problem_best_parameter_print` */
static coco_problem_stacked_t *coco_problem_as_stacked(const coco_problem_t *problem);

void bbob_biobj_problem_best_parameter_print(const coco_problem_t *problem) {
  coco_file_t file;
  coco_problem_t *problem1, *problem2;

  /* We need the underlying single objective problems of the stacked problem, which are not accessible
   * once the problem has been wrapped, e.g. with a logger.
   */
  if (dynamic_cast<const coco_problem_stacked_t *>(problem) == nullptr) {
    coco_warning("I'm sorry, Dave. I'm afraid I can't do that.");
    /* In case this function was called previously, overwrite any results. */
    file.open("._bbob_biobj_problem_best_parameter.txt", "w");
    return;
  }

  assert(problem != nullptr);
  problem1 = coco_problem_as_stacked(problem)->get_problem1();
  assert(problem1 != nullptr);
  assert(!problem1->get_best_parameter().empty());

  problem2 = coco_problem_as_stacked(problem)->get_problem2();
  assert(problem2 != nullptr);
  assert(!problem2->get_best_parameter().empty());

  if (file.open("._bbob_biobj_problem_best_parameter.txt", "w")) {
    for (size_t i = 0; i < problem->get_number_of_variables(); ++i)
      fprintf(file, " %.16f ", problem1->get_best_parameter()[i]);
    fprintf(file, "\n");
    for (size_t i = 0; i < problem->get_number_of_variables(); ++i)
      fprintf(file, " %.16f ", problem2->get_best_parameter()[i]);
  }
}

/***********************************************************************************************************/

/**
 * @name Methods regarding the transformed COCO problem
 */
/**@{*/

/**
 * @brief Returns the transformed problem behind the given problem pointer.
 */
static coco_problem_transformed_t *coco_problem_as_transformed(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(dynamic_cast<const coco_problem_transformed_t *>(problem) != nullptr);
  return static_cast<coco_problem_transformed_t *>(const_cast<coco_problem_t *>(problem));
}

/**
 * @brief Constructs a transformed problem that wraps the inner_problem.
 *
 * By default all methods will dispatch to the inner_problem. A prefix is prepended to the problem name
 * in order to reflect the transformation somewhere.
 */
coco_problem_transformed_s::coco_problem_transformed_s(coco_problem_ptr_t inner, std::string_view name_prefix)
    : coco_problem_s(inner->get_number_of_variables(), inner->get_number_of_objectives(),
                     inner->get_number_of_constraints()),
      inner_problem(std::move(inner)) {
  copy_properties_from(*inner_problem);
  set_name("%.*s(%s)", static_cast<int>(name_prefix.size()), name_prefix.data(), inner_problem->get_name().c_str());
}

/**
 * @brief Evaluates the objective function of the inner problem.
 */
void coco_problem_transformed_s::evaluate_function(std::span<const double> x, std::span<double> y) {
  inner_problem->evaluate_function(x, y);
}

/**
 * @brief Evaluates the constraints of the inner problem.
 */
void coco_problem_transformed_s::evaluate_constraint(std::span<const double> x, std::span<double> y,
                                                     int update_counter) {
  inner_problem->evaluate_constraint(x, y, update_counter);
}

/**
 * @brief Calls the bbob_evaluate_gradient function on the inner problem.
 */
void coco_problem_transformed_s::evaluate_gradient(std::span<const double> x, std::span<double> y) {
  inner_problem->evaluate_gradient(x, y);
}

/**
 * @brief Calls the coco_recommend_solution function on the inner problem.
 */
void coco_problem_transformed_s::recommend_solution(std::span<const double> x) { inner_problem->recommend_solution(x); }

/**
 * @brief Maps base_x to the coordinates of the inner problem and from there to the coordinates of this problem.
 */
bool coco_problem_transformed_s::from_base_coordinates(std::span<const double> base_x, std::span<double> x) const {
  std::vector<double> inner_x(inner_problem->get_number_of_variables());
  return inner_problem->from_base_coordinates(base_x, inner_x) && inverse(inner_x, x);
}

void coco_problem_transformed_s::update_best_parameter() {
  if (best_parameter.empty())
    return;
  if (!inverse(inner_problem->get_best_parameter(), best_parameter)) {
    coco_warning("%s: the transformation is not invertible, 'best_parameter' set to NAN", get_name().c_str());
    coco_vector_set_to_nan(best_parameter);
  }
}

bool coco_problem_transformed_s::identity_inverse(std::span<const double> inner_x, std::span<double> x) const {
  assert(inner_x.size() == x.size());
  std::copy(inner_x.begin(), inner_x.end(), x.begin());
  return true;
}

/**
 * @brief Frees the transformed problem but not its inner problem, which is returned instead.
 */
static coco_problem_ptr_t coco_problem_transformed_free_keep_inner(coco_problem_ptr_t problem) {
  return coco_problem_as_transformed(problem.get())->release_inner_problem();
}

/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods regarding the stacked COCO problem
 */
/**@{*/

/**
 * @brief Returns the stacked problem behind the given problem pointer.
 */
static coco_problem_stacked_t *coco_problem_as_stacked(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(dynamic_cast<const coco_problem_stacked_t *>(problem) != nullptr);
  return static_cast<coco_problem_stacked_t *>(const_cast<coco_problem_t *>(problem));
}

coco_problem_stacked_s::coco_problem_stacked_s(coco_problem_ptr_t first, coco_problem_ptr_t second,
                                               std::span<const double> smallest_values_of_interest,
                                               std::span<const double> largest_values_of_interest)
    : coco_problem_s(
          coco_problem_get_dimension(first.get()),
          coco_problem_get_number_of_objectives(first.get()) + coco_problem_get_number_of_objectives(second.get()),
          coco_problem_get_number_of_constraints(first.get()) + coco_problem_get_number_of_constraints(second.get())),
      problem1(std::move(first)), problem2(std::move(second)) {

  assert(problem1);
  assert(problem2);
  assert(coco_problem_get_dimension(problem1.get()) == coco_problem_get_dimension(problem2.get()));

  set_id("%s__%s", problem1->get_id().c_str(), problem2->get_id().c_str());
  set_name("%s + %s", problem1->get_name().c_str(), problem2->get_name().c_str());

  /* Without an explicit region of interest, use the one of the second problem */
  if (smallest_values_of_interest.empty())
    smallest_values_of_interest = problem2->get_smallest_values_of_interest();
  if (largest_values_of_interest.empty())
    largest_values_of_interest = problem2->get_largest_values_of_interest();

  for (size_t i = 0; i < number_of_variables; ++i) {
    this->smallest_values_of_interest[i] = smallest_values_of_interest[i];
    this->largest_values_of_interest[i] = largest_values_of_interest[i];
  }
  assert(problem1->get_number_of_integer_variables() == problem2->get_number_of_integer_variables());
  number_of_integer_variables = problem1->get_number_of_integer_variables();

  assert(!best_value.empty());

  if (number_of_constraints > 0) {

    /* The best_value must be set up afterwards in suite_cons_bbob_problems.cpp */
    best_value[0] = -FLT_MAX;

    /* Define initial_solution as problem2->initial_solution */
    if (coco_problem_get_number_of_constraints(problem2.get()) > 0 && !problem2->get_initial_solution().empty())
      initial_solution = problem2->get_initial_solution();

  } else {

    /* Compute the ideal and nadir values */
    assert(!nadir_value.empty());

    best_value[0] = problem1->get_best_value()[0];
    best_value[1] = problem2->get_best_value()[0];
    coco_evaluate_function(problem1.get(), problem2->get_best_parameter().data(), &nadir_value[0]);
    coco_evaluate_function(problem2.get(), problem1->get_best_parameter().data(), &nadir_value[1]);
  }
}

/**
 * @brief Evaluates the objective functions of the underlying problems.
 */
void coco_problem_stacked_s::evaluate_function(std::span<const double> x, std::span<double> y) {
  const size_t number_of_objectives_problem1 = coco_problem_get_number_of_objectives(problem1.get());
  const size_t number_of_objectives_problem2 = coco_problem_get_number_of_objectives(problem2.get());
  assert(coco_problem_get_number_of_objectives(this) == number_of_objectives_problem1 + number_of_objectives_problem2);

  if (number_of_objectives_problem1 > 0)
    problem1->evaluate_function(x, y.first(number_of_objectives_problem1));
  if (number_of_objectives_problem2 > 0)
    problem2->evaluate_function(x, y.subspan(number_of_objectives_problem1));

  /* Make sure that no feasible point has a function value lower
   * than the minimum's.
   */
  if (number_of_constraints > 0)
    assert_not_below_best_value(x, y);
}

/**
 * @brief Evaluates the constraints of the underlying problems.
 */
void coco_problem_stacked_s::evaluate_constraint(std::span<const double> x, std::span<double> y, int update_counter) {
  const size_t number_of_constraints_problem1 = coco_problem_get_number_of_constraints(problem1.get());
  const size_t number_of_constraints_problem2 = coco_problem_get_number_of_constraints(problem2.get());

  if (number_of_constraints == 0) {
    coco_problem_s::evaluate_constraint(x, y, update_counter);
    return;
  }

  assert(coco_problem_get_number_of_constraints(this) ==
         number_of_constraints_problem1 + number_of_constraints_problem2);

  if (number_of_constraints_problem1 > 0)
    problem1->evaluate_constraint(x, y.first(number_of_constraints_problem1), update_counter);
  if (number_of_constraints_problem2 > 0)
    problem2->evaluate_constraint(x, y.subspan(number_of_constraints_problem1), update_counter);
}

/* TODO: Missing coco_problem_stacked_recommend_solution function! */

/**
 * @brief Allocates a problem constructed by stacking two COCO problems.
 *
 * This is particularly useful for generating multi-objective problems, e.g. a bi-objective problem from two
 * single-objective problems. The stacked problem must behave like a normal COCO problem accepting the same
 * input. Without a region of interest, the one of the second problem is used.
 *
 * @note Best parameter becomes somewhat meaningless, but the nadir value make sense now.
 */
static coco_problem_ptr_t coco_problem_stacked_allocate(coco_problem_ptr_t first, coco_problem_ptr_t second,
                                                        std::span<const double> smallest_values_of_interest = {},
                                                        std::span<const double> largest_values_of_interest = {}) {
  return std::make_unique<coco_problem_stacked_t>(std::move(first), std::move(second), smallest_values_of_interest,
                                                  largest_values_of_interest);
}
/**@}*/

/***********************************************************************************************************/
