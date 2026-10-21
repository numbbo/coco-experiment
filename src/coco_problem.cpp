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
  size_t i, j;
  int is_feasible;

  assert(problem != nullptr);

  /* Set objective vector to INFINITY if the decision vector contains any INFINITY values */
  for (i = 0; i < coco_problem_get_dimension(problem); i++) {
    if (coco_is_inf(x[i])) {
      for (j = 0; j < coco_problem_get_number_of_objectives(problem); j++) {
        y[j] = fabs(x[i]);
      }
      return;
    }
  }

  /* Set objective vector to NAN if the decision vector contains any NAN values */
  if (coco_vector_contains_nan(x, coco_problem_get_dimension(problem))) {
    coco_vector_set_to_nan(y, coco_problem_get_number_of_objectives(problem));
    return;
  }

  problem->evaluate_function(x, y);
  problem->evaluations++; /* each derived class has its own counter, only the most outer will be visible */

  /* A little bit of bookkeeping */
  if (y[0] < problem->best_observed_fvalue[0]) {
    is_feasible = 1;
    if (coco_problem_get_number_of_constraints(problem) > 0) {
      std::vector<double> z(coco_problem_get_number_of_constraints(problem));
      is_feasible = coco_is_feasible(problem, x, z.data());
    }
    if (is_feasible) {
      problem->best_observed_fvalue[0] = y[0];
      problem->best_observed_evaluation[0] = problem->evaluations;
    }
  }
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
  size_t i, j;
  assert(problem != nullptr);

  /* Set constraints vector to INFINITY if the decision vector contains any INFINITY values */
  for (i = 0; i < coco_problem_get_dimension(problem); i++) {
    if (coco_is_inf(x[i])) {
      for (j = 0; j < coco_problem_get_number_of_constraints(problem); j++) {
        y[j] = fabs(x[i]);
      }
      return;
    }
  }

  /* Set constraints vector to NAN if the decision vector contains any NAN values */
  if (coco_vector_contains_nan(x, coco_problem_get_dimension(problem))) {
    coco_vector_set_to_nan(y, coco_problem_get_number_of_constraints(problem));
    return;
  }

  problem->evaluate_constraint(x, y, update_counter);
  if (update_counter)
    problem->evaluations_constraints++;
}

void coco_evaluate_constraint(coco_problem_t *problem, const double *x, double *y) {
  coco_evaluate_constraint_optional_update(problem, x, y, 1);
}

/**
 * @note Both x and y must point to correctly sized allocated memory regions.
 *
 * @param problem The given COCO problem.
 * @param x The decision vector.
 * @param y The gradient of the function evaluated at the point x.
 */
static void bbob_evaluate_gradient(coco_problem_t *problem, const double *x, double *y) {
  assert(problem != nullptr);
  problem->evaluate_gradient(x, y);
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
  problem->recommend_solution(x);
}

/**
 * @name Default implementations of the problem methods
 *
 * They report that the method is not implemented for the problem.
 */
/**@{*/
void coco_problem_s::evaluate_function(const double *x, double *y) {
  (void)x; /* To silence the compiler */
  (void)y;
  coco_error("coco_evaluate_function(): No function implemented for problem %s", problem_id.c_str());
}

void coco_problem_s::evaluate_constraint(const double *x, double *y, int update_counter) {
  (void)x; /* To silence the compiler */
  (void)y;
  (void)update_counter;
  coco_error("coco_evaluate_constraint_optional_update(): No constraint function implemented for problem %s",
             problem_id.c_str());
}

void coco_problem_s::evaluate_gradient(const double *x, double *y) {
  (void)x; /* To silence the compiler */
  (void)y;
  coco_error("bbob_evaluate_gradient(): No gradient function implemented for problem %s", problem_id.c_str());
}

void coco_problem_s::recommend_solution(const double *x) {
  (void)x; /* To silence the compiler */
  coco_warning("coco_recommend_solution(): No recommend solution function implemented for problem %s",
               problem_id.c_str());
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

  is_noisy = 0;
  if (number_of_objectives > 0)
    last_noise_free_values.assign(number_of_objectives, 0.0);
  if (number_of_objectives > 1) {
    is_opt_known = 0; /* Optimum of multi-objective problems is unknown by default */
    best_value.assign(number_of_objectives, 0.0);
    nadir_value.assign(number_of_objectives, 0.0);
  } else {
    is_opt_known = 1; /* Optimum of single-objective problems is known by default */
    best_parameter.assign(number_of_variables, 0.0);
    best_value.assign(1, 0.0);
  }
  is_tainted = 0;

  evaluations = 0;
  evaluations_constraints = 0;
  final_target_delta[0] = 1e-8; /* in case to be modified by the benchmark */
  best_observed_fvalue[0] = DBL_MAX;
  best_observed_evaluation[0] = 0;
  suite = nullptr; /* To be initialized in the coco_suite_get_problem_from_indices() function */
  suite_dep_index = 0;
  suite_dep_function = 0;
  suite_dep_instance = 0;
}

/**
 * @brief Allocates a new coco_problem_t for the given number of variables, number of objectives and
 * number of constraints.
 */
static coco_problem_ptr_t coco_problem_allocate(const size_t number_of_variables, const size_t number_of_objectives,
                                                const size_t number_of_constraints) {
  return std::make_unique<coco_problem_t>(number_of_variables, number_of_objectives, number_of_constraints);
}

/**
 * @brief Copies the properties of the 'other' problem, i.e. everything except its behaviour, data and noise.
 *
 * Both problems must have the same number of variables, objectives and constraints.
 */
static void coco_problem_copy_properties(coco_problem_t *problem, const coco_problem_t *other) {
  problem->smallest_values_of_interest = other->smallest_values_of_interest;
  problem->largest_values_of_interest = other->largest_values_of_interest;
  if (!other->best_parameter.empty())
    problem->best_parameter = other->best_parameter;
  problem->is_tainted = other->is_tainted;

  problem->number_of_integer_variables = other->number_of_integer_variables;

  problem->initial_solution = other->initial_solution;

  problem->is_opt_known = other->is_opt_known;
  if (!other->best_value.empty())
    problem->best_value = other->best_value;
  if (!other->nadir_value.empty())
    problem->nadir_value = other->nadir_value;

  problem->problem_name = other->problem_name;
  problem->problem_id = other->problem_id;
  problem->problem_type = other->problem_type;

  problem->evaluations = other->evaluations;
  problem->evaluations_constraints = other->evaluations_constraints;
  problem->final_target_delta[0] = other->final_target_delta[0];
  problem->best_observed_fvalue[0] = other->best_observed_fvalue[0];
  problem->best_observed_evaluation[0] = other->best_observed_evaluation[0];

  problem->suite = other->suite;
  problem->suite_dep_index = other->suite_dep_index;
  problem->suite_dep_function = other->suite_dep_function;
  problem->suite_dep_instance = other->suite_dep_instance;
}

/**
 * @brief Initializes a single-objective problem using scalar values for smallest_value_of_interest,
 * largest_value_of_interest and best_parameter. Assumes all variables are continuous.
 */
static void coco_problem_init_from_scalars(coco_problem_t *problem, std::string_view problem_name,
                                           const double smallest_value_of_interest,
                                           const double largest_value_of_interest, const double best_parameter) {
  size_t i;

  assert(problem->number_of_objectives == 1);
  assert(problem->number_of_constraints == 0);
  problem->problem_name = problem_name;

  for (i = 0; i < problem->number_of_variables; ++i) {
    problem->smallest_values_of_interest[i] = smallest_value_of_interest;
    problem->largest_values_of_interest[i] = largest_value_of_interest;
    problem->best_parameter[i] = best_parameter;
  }
  problem->number_of_integer_variables = 0;
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
static int coco_problem_id_is_fine(const std::string &id) {
  for (const char c : id) {
    if (('A' <= c) && (c <= 'Z'))
      continue;
    if (('a' <= c) && (c <= 'z'))
      continue;
    if ((c == '_') || (c == '-'))
      continue;
    if (('0' <= c) && (c <= '9'))
      continue;
    return 0;
  }
  return 1;
}

/**
 * @brief Sets the problem_id using formatted printing (as in printf).
 *
 * Verifies that the problem_id is in the correct format.
 */
COCO_FORMAT(2, 3) static void coco_problem_set_id(coco_problem_t *problem, const char *id, ...) {
  va_list args;

  va_start(args, id);
  problem->problem_id = coco_string_vformat(id, args);
  va_end(args);
  if (!coco_problem_id_is_fine(problem->problem_id)) {
    coco_error("Problem id should only contain standard chars, not like '%s'", problem->problem_id.c_str());
  }
}

/**
 * @brief Sets the problem_name using formatted printing (as in printf).
 */
COCO_FORMAT(2, 3) static void coco_problem_set_name(coco_problem_t *problem, const char *name, ...) {
  va_list args;

  va_start(args, name);
  problem->problem_name = coco_string_vformat(name, args);
  va_end(args);
}

/**
 * @brief Sets the problem_type using formatted printing (as in printf).
 */
COCO_FORMAT(2, 3) static void coco_problem_set_type(coco_problem_t *problem, const char *type, ...) {
  va_list args;

  va_start(args, type);
  problem->problem_type = coco_string_vformat(type, args);
  va_end(args);
}

size_t coco_problem_get_evaluations(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return problem->evaluations;
}

size_t coco_problem_get_evaluations_constraints(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return problem->evaluations_constraints;
}

/**
 * @brief Returns 1 if the best parameter is not (close to) zero and 0 otherwise.
 */
static int coco_problem_best_parameter_not_zero(const coco_problem_t *problem) {
  size_t i = 0;
  int best_is_zero = 1;

  if (coco_vector_contains_nan(problem->best_parameter.data(), problem->number_of_variables))
    return 1;

  while (i < problem->number_of_variables && best_is_zero) {
    best_is_zero = coco_double_almost_equal(problem->best_parameter[i], 0, 1e-9);
    i++;
  }

  return !best_is_zero;
}

/**
 * @note Can be used to prevent unnecessary burning of CPU time.
 */
int coco_problem_final_target_hit(const coco_problem_t *problem) {
  assert(problem != nullptr);
  if (coco_problem_get_number_of_objectives(problem) != 1 || coco_problem_get_evaluations(problem) < 1)
    return 0;
  if (problem->best_value.empty())
    return 0;
  return problem->best_observed_fvalue[0] <= problem->best_value[0] + problem->final_target_delta[0] ? 1 : 0;
}
/**
 * @note Tentative...
 */
double coco_problem_get_best_observed_fvalue1(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return problem->best_observed_fvalue[0];
}

/**
 * @brief Returns the optimal function value of the problem
 */
double coco_problem_get_best_value(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(!problem->best_value.empty());
  return problem->best_value[0];
}

/**
 * @note This function breaks the black-box property: the returned  value is not
 * meant to be used by the optimization algorithm.
 */
double coco_problem_get_final_target_fvalue1(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(!problem->best_value.empty());
  return problem->best_value[0] + problem->final_target_delta[0];
}

/**
 * @note Do not modify the returned string! If you free the problem, the returned pointer becomes invalid.
 * When in doubt, copy the returned value.
 */
const char *coco_problem_get_name(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(!problem->problem_name.empty());
  return problem->problem_name.c_str();
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
  assert(!problem->problem_id.empty());
  return problem->problem_id.c_str();
}

const char *coco_problem_get_type(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(!problem->problem_type.empty());
  return problem->problem_type.c_str();
}

size_t coco_problem_get_dimension(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(problem->number_of_variables > 0);
  return problem->number_of_variables;
}

size_t coco_problem_get_number_of_objectives(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return problem->number_of_objectives;
}

size_t coco_problem_get_number_of_constraints(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return problem->number_of_constraints;
}

const double *coco_problem_get_smallest_values_of_interest(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(!problem->smallest_values_of_interest.empty());
  return problem->smallest_values_of_interest.data();
}

const double *coco_problem_get_largest_values_of_interest(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(!problem->largest_values_of_interest.empty());
  return problem->largest_values_of_interest.data();
}

size_t coco_problem_get_number_of_integer_variables(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return problem->number_of_integer_variables;
}

const double *coco_problem_get_largest_fvalues_of_interest(const coco_problem_t *problem) {
  assert(problem != nullptr);
  if (problem->number_of_objectives == 1)
    coco_error(
        "coco_problem_get_largest_fvalues_of_interest(): f-values of interest undefined for single-objective problems");
  if (problem->nadir_value.empty())
    coco_error("coco_problem_get_largest_fvalues_of_interest(): f-values of interest undefined");
  return problem->nadir_value.data();
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

  size_t i;

  assert(problem != nullptr);
  if (!problem->initial_solution.empty()) {
    for (i = 0; i < problem->number_of_variables; ++i)
      initial_solution[i] = problem->initial_solution[i];
  } else {
    assert(!problem->smallest_values_of_interest.empty());
    assert(!problem->largest_values_of_interest.empty());
    for (i = 0; i < problem->number_of_variables; ++i)
      initial_solution[i] = problem->smallest_values_of_interest[i] +
                            0.5 * (problem->largest_values_of_interest[i] - problem->smallest_values_of_interest[i]);
    if (problem->number_of_integer_variables > 0) {
      for (i = 0; i < problem->number_of_integer_variables; ++i) {
        initial_solution[i] = coco_double_round(initial_solution[i]);
      }
    }
  }
}

static coco_suite_t *coco_problem_get_suite(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return problem->suite;
}

static void coco_problem_set_suite(coco_problem_t *problem, coco_suite_t *suite) {
  assert(problem != nullptr);
  problem->suite = suite;
}

size_t coco_problem_get_suite_dep_index(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return problem->suite_dep_index;
}

static size_t coco_problem_get_suite_dep_function(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(problem->suite_dep_function > 0);
  return problem->suite_dep_function;
}

static size_t coco_problem_get_suite_dep_instance(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(problem->suite_dep_instance > 0);
  return problem->suite_dep_instance;
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
  size_t i;
  assert(problem != nullptr);

  if (problem->best_parameter.empty() || problem->number_of_objectives != 1) {
    return 0;
  }

  // No more benchmarking on this problem!
  problem->is_tainted = 1;

  for (i = 0; i < problem->number_of_variables; ++i) {
    best_parameter[i] = problem->best_parameter[i];
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
int coco_problem_is_tainted(const coco_problem_t *problem) { return problem->is_tainted; }

void bbob_problem_best_parameter_print(const coco_problem_t *problem) {
  size_t i;
  FILE *file;
  assert(problem != nullptr);
  assert(!problem->best_parameter.empty());
  file = fopen("._bbob_problem_best_parameter.txt", "w");
  if (file != nullptr) {
    for (i = 0; i < problem->number_of_variables; ++i)
      fprintf(file, " %.16f ", problem->best_parameter[i]);
    fclose(file);
  }
}

/* We need a forward declaration of this for `bbob_biobj_problem_best_parameter_print` */
static coco_problem_stacked_t *coco_problem_as_stacked(const coco_problem_t *problem);

void bbob_biobj_problem_best_parameter_print(const coco_problem_t *problem) {
  size_t i;
  FILE *file;
  coco_problem_t *problem1, *problem2;

  /* We need the underlying single objective problems of the stacked problem, which are not accessible
   * once the problem has been wrapped, e.g. with a logger.
   */
  if (dynamic_cast<const coco_problem_stacked_t *>(problem) == nullptr) {
    coco_warning("I'm sorry, Dave. I'm afraid I can't do that.");
    /* In case this function was called previously, overwrite any results. */
    file = fopen("._bbob_biobj_problem_best_parameter.txt", "w");
    if (file != nullptr) {
      fclose(file);
    }
    return;
  }

  assert(problem != nullptr);
  problem1 = coco_problem_as_stacked(problem)->problem1.get();
  assert(problem1 != nullptr);
  assert(!problem1->best_parameter.empty());

  problem2 = coco_problem_as_stacked(problem)->problem2.get();
  assert(problem2 != nullptr);
  assert(!problem2->best_parameter.empty());

  file = fopen("._bbob_biobj_problem_best_parameter.txt", "w");
  if (file != nullptr) {
    for (i = 0; i < problem->number_of_variables; ++i)
      fprintf(file, " %.16f ", problem1->best_parameter[i]);
    fprintf(file, "\n");
    for (i = 0; i < problem->number_of_variables; ++i)
      fprintf(file, " %.16f ", problem2->best_parameter[i]);
    fclose(file);
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
 * @brief Returns the inner problem of the transformed problem.
 */
static coco_problem_t *coco_problem_transformed_get_inner_problem(const coco_problem_t *problem) {
  coco_problem_transformed_t *transformed = coco_problem_as_transformed(problem);
  assert(transformed->inner_problem != nullptr);
  return transformed->inner_problem.get();
}

/**
 * @brief Constructs a transformed problem that wraps the inner_problem.
 *
 * By default all methods will dispatch to the inner_problem. A prefix is prepended to the problem name
 * in order to reflect the transformation somewhere.
 */
coco_problem_transformed_s::coco_problem_transformed_s(coco_problem_ptr_t inner, std::string_view name_prefix)
    : coco_problem_s(inner->number_of_variables, inner->number_of_objectives, inner->number_of_constraints),
      inner_problem(std::move(inner)) {
  coco_problem_copy_properties(this, inner_problem.get());
  problem_name = std::string(name_prefix) + "(" + inner_problem->problem_name + ")";
}

/**
 * @brief Calls the coco_evaluate_function function on the inner problem.
 */
void coco_problem_transformed_s::evaluate_function(const double *x, double *y) {
  coco_evaluate_function(inner_problem.get(), x, y);
}

/**
 * @brief Calls the coco_evaluate_constraint_optional_update function on the inner problem.
 */
void coco_problem_transformed_s::evaluate_constraint(const double *x, double *y, int update_counter) {
  coco_evaluate_constraint_optional_update(inner_problem.get(), x, y, update_counter);
}

/**
 * @brief Calls the bbob_evaluate_gradient function on the inner problem.
 */
void coco_problem_transformed_s::evaluate_gradient(const double *x, double *y) {
  bbob_evaluate_gradient(inner_problem.get(), x, y);
}

/**
 * @brief Calls the coco_recommend_solution function on the inner problem.
 */
void coco_problem_transformed_s::recommend_solution(const double *x) { coco_recommend_solution(inner_problem.get(), x); }

/**
 * @brief Frees the transformed problem but not its inner problem, which is returned instead.
 */
static coco_problem_ptr_t coco_problem_transformed_free_keep_inner(coco_problem_ptr_t problem) {
  return std::move(coco_problem_as_transformed(problem.get())->inner_problem);
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

coco_problem_stacked_s::coco_problem_stacked_s(coco_problem_ptr_t first, coco_problem_ptr_t second)
    : coco_problem_s(
          coco_problem_get_dimension(first.get()),
          coco_problem_get_number_of_objectives(first.get()) + coco_problem_get_number_of_objectives(second.get()),
          coco_problem_get_number_of_constraints(first.get()) + coco_problem_get_number_of_constraints(second.get())),
      problem1(std::move(first)), problem2(std::move(second)) {}

/**
 * @brief Calls the coco_evaluate_function function on the underlying problems.
 */
void coco_problem_stacked_s::evaluate_function(const double *x, double *y) {
  const size_t number_of_objectives_problem1 = coco_problem_get_number_of_objectives(problem1.get());
  const size_t number_of_objectives_problem2 = coco_problem_get_number_of_objectives(problem2.get());
  int is_feasible;

  assert(coco_problem_get_number_of_objectives(this) == number_of_objectives_problem1 + number_of_objectives_problem2);

  if (number_of_objectives_problem1 > 0)
    coco_evaluate_function(problem1.get(), x, &y[0]);
  if (number_of_objectives_problem2 > 0)
    coco_evaluate_function(problem2.get(), x, &y[number_of_objectives_problem1]);

  /* Make sure that no feasible point has a function value lower
   * than the minimum's.
   */
  if (number_of_constraints > 0) {
    std::vector<double> cons_values(number_of_constraints);
    is_feasible = coco_is_feasible(this, x, cons_values.data());
    if (is_feasible)
      assert(y[0] + 1e-13 >= best_value[0]);
  }
}

/**
 * @brief Calls the coco_evaluate_constraint_optional_update function on the underlying problems.
 */
void coco_problem_stacked_s::evaluate_constraint(const double *x, double *y, int update_counter) {
  const size_t number_of_constraints_problem1 = coco_problem_get_number_of_constraints(problem1.get());
  const size_t number_of_constraints_problem2 = coco_problem_get_number_of_constraints(problem2.get());

  if (number_of_constraints == 0) {
    coco_problem_s::evaluate_constraint(x, y, update_counter);
    return;
  }

  assert(coco_problem_get_number_of_constraints(this) ==
         number_of_constraints_problem1 + number_of_constraints_problem2);

  if (number_of_constraints_problem1 > 0)
    coco_evaluate_constraint_optional_update(problem1.get(), x, y, update_counter);
  if (number_of_constraints_problem2 > 0)
    coco_evaluate_constraint_optional_update(problem2.get(), x, &y[number_of_constraints_problem1], update_counter);
}

/* TODO: Missing coco_problem_stacked_recommend_solution function! */

/**
 * @brief Allocates a problem constructed by stacking two COCO problems.
 *
 * This is particularly useful for generating multi-objective problems, e.g. a bi-objective problem from two
 * single-objective problems. The stacked problem must behave like a normal COCO problem accepting the same
 * input.
 *
 * @note Regions of interest in the decision space must either agree or at least one of them must be NULL.
 * @note Best parameter becomes somewhat meaningless, but the nadir value make sense now.
 */
static coco_problem_ptr_t coco_problem_stacked_allocate(coco_problem_ptr_t first, coco_problem_ptr_t second,
                                                        const double *smallest_values_of_interest,
                                                        const double *largest_values_of_interest) {

  size_t number_of_variables, number_of_constraints;
  size_t i;

  assert(first);
  assert(second);
  assert(coco_problem_get_dimension(first.get()) == coco_problem_get_dimension(second.get()));

  auto stacked = std::make_unique<coco_problem_stacked_t>(std::move(first), std::move(second));
  coco_problem_t *problem = stacked.get(); /* the new coco problem */
  const coco_problem_t *problem1 = stacked->problem1.get();
  const coco_problem_t *problem2 = stacked->problem2.get();
  number_of_variables = problem->number_of_variables;
  number_of_constraints = problem->number_of_constraints;

  problem->problem_id = problem1->problem_id + "__" + problem2->problem_id;
  problem->problem_name = problem1->problem_name + " + " + problem2->problem_name;

  assert(smallest_values_of_interest);
  assert(largest_values_of_interest);

  for (i = 0; i < number_of_variables; ++i) {
    problem->smallest_values_of_interest[i] = smallest_values_of_interest[i];
    problem->largest_values_of_interest[i] = largest_values_of_interest[i];
  }
  assert(problem1->number_of_integer_variables == problem2->number_of_integer_variables);
  problem->number_of_integer_variables = problem1->number_of_integer_variables;

  assert(!problem->best_value.empty());

  if (number_of_constraints > 0) {

    /* The best_value must be set up afterwards in suite_cons_bbob_problems.cpp */
    problem->best_value[0] = -FLT_MAX;

    /* Define problem->initial_solution as problem2->initial_solution */
    if (coco_problem_get_number_of_constraints(problem2) > 0 && !problem2->initial_solution.empty())
      problem->initial_solution = problem2->initial_solution;

  } else {

    /* Compute the ideal and nadir values */
    assert(!problem->nadir_value.empty());

    problem->best_value[0] = problem1->best_value[0];
    problem->best_value[1] = problem2->best_value[0];
    coco_evaluate_function(stacked->problem1.get(), problem2->best_parameter.data(), &problem->nadir_value[0]);
    coco_evaluate_function(stacked->problem2.get(), problem1->best_parameter.data(), &problem->nadir_value[1]);
  }

  return stacked;
}

/**
 * @brief Allocates a problem constructed by stacking two COCO problems, using the region of interest of the
 * second problem.
 */
static coco_problem_ptr_t coco_problem_stacked_allocate(coco_problem_ptr_t first, coco_problem_ptr_t second) {
  const double *smallest_values_of_interest = second->smallest_values_of_interest.data();
  const double *largest_values_of_interest = second->largest_values_of_interest.data();
  return coco_problem_stacked_allocate(std::move(first), std::move(second), smallest_values_of_interest,
                                       largest_values_of_interest);
}
/**@}*/

/***********************************************************************************************************/
