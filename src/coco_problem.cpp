/**
 * @file coco_problem.cpp
 * @brief Definitions of functions regarding COCO problems.
 */

#include <float.h>
#include "coco.h"
#include "coco_internal.h"

#include "coco_utilities.cpp"

namespace coco {

/***********************************************************************************************************/

/**
 * @name Methods regarding the basic COCO problem
 */
/**@{*/
bool Problem::SetNonFiniteValues(std::span<const double> x, std::span<double> y) {
  /* Set y to INFINITY if x contains any INFINITY values */
  for (const double xi : x) {
    if (IsInf(xi)) {
      for (double &yj : y) {
        yj = fabs(xi);
      }
      return true;
    }
  }

  /* Set y to NAN if x contains any NAN values */
  if (VectorContainsNan(x)) {
    VectorSetToNan(y);
    return true;
  }
  return false;
}

void Problem::EvaluateFunctionChecked(std::span<const double> x, std::span<double> y) {
  if (SetNonFiniteValues(x, y))
    return;
  EvaluateFunction(x, y);
  CountEvaluation(x, y);
}

void Problem::EvaluateConstraintChecked(std::span<const double> x, std::span<double> y, int update_counter) {
  if (SetNonFiniteValues(x, y))
    return;
  EvaluateConstraint(x, y, update_counter);
  if (update_counter)
    CountConstraintEvaluation();
}

/**
 * @name Default implementations of the problem methods
 *
 * They report that the method is not implemented for the problem.
 */
/**@{*/
void Problem::EvaluateFunction(std::span<const double>, std::span<double>) {
  coco_error("coco_evaluate_function(): No function implemented for problem %s", problem_id_.c_str());
}

void Problem::EvaluateConstraint(std::span<const double>, std::span<double>, int) {
  coco_error("EvaluateConstraint(): No constraint function implemented for problem %s", problem_id_.c_str());
}

void Problem::EvaluateGradient(std::span<const double>, std::span<double>) {
  coco_error("EvaluateGradient(): No gradient function implemented for problem %s", problem_id_.c_str());
}

void Problem::RecommendSolution(std::span<const double>) {
  coco_warning("coco_recommend_solution(): No recommend solution function implemented for problem %s",
               problem_id_.c_str());
}

void Problem::AssertNotBelowBestValue([[maybe_unused]] std::span<const double> x,
                                      [[maybe_unused]] std::span<const double> y) {
#ifndef NDEBUG
  if (number_of_constraints_ == 0 || IsFeasible(x, {}))
    assert(y[0] + 1e-13 >= best_value_[0]);
#endif
}

bool Problem::FromBaseCoordinates(std::span<const double> base_x, std::span<double> x) const {
  assert(base_x.size() == x.size());
  std::copy(base_x.begin(), base_x.end(), x.begin());
  return true;
}
/**@}*/
/***********************************************************************************************************/

/**
 * @brief Constructs a problem for the given number of variables, number of objectives and number of constraints.
 */
Problem::Problem(const size_t number_of_variables, const size_t number_of_objectives,
                 const size_t number_of_constraints) {
  /* Initialize fields to sane/safe defaults */
  this->number_of_variables_ = number_of_variables;
  this->number_of_objectives_ = number_of_objectives;
  this->number_of_constraints_ = number_of_constraints;
  smallest_values_of_interest_.assign(number_of_variables, 0.0);
  largest_values_of_interest_.assign(number_of_variables, 0.0);
  number_of_integer_variables_ = 0; /* No integer variables by default */

  noisy_ = false;
  if (number_of_objectives > 0)
    last_noise_free_values_.assign(number_of_objectives, 0.0);
  if (number_of_objectives > 1) {
    opt_known_ = false; /* Optimum of multi-objective problems is unknown by default */
    best_value_.assign(number_of_objectives, 0.0);
    nadir_value_.assign(number_of_objectives, 0.0);
  } else {
    opt_known_ = true; /* Optimum of single-objective problems is known by default */
    best_parameter_.assign(number_of_variables, 0.0);
    best_value_.assign(1, 0.0);
  }
  tainted_ = false;

  evaluations_ = 0;
  evaluations_constraints_ = 0;
  final_target_delta_ = 1e-8; /* in case to be modified by the benchmark */
  best_observed_fvalue_ = DBL_MAX;
  best_observed_evaluation_ = 0;
  suite_ = nullptr; /* To be initialized in the coco_suite_get_problem_from_indices() function */
  suite_dep_index_ = 0;
  suite_dep_function_ = 0;
  suite_dep_instance_ = 0;
}

/**
 * @brief Copies the properties of the 'other' problem, i.e. everything except its behaviour, data and noise.
 *
 * Both problems must have the same number of variables, objectives and constraints.
 */
void Problem::CopyPropertiesFrom(const Problem &other) {
  smallest_values_of_interest_ = other.smallest_values_of_interest_;
  largest_values_of_interest_ = other.largest_values_of_interest_;
  if (!other.best_parameter_.empty())
    best_parameter_ = other.best_parameter_;
  tainted_ = other.tainted_;

  number_of_integer_variables_ = other.number_of_integer_variables_;

  initial_solution_ = other.initial_solution_;

  opt_known_ = other.opt_known_;
  if (!other.best_value_.empty())
    best_value_ = other.best_value_;
  if (!other.nadir_value_.empty())
    nadir_value_ = other.nadir_value_;

  problem_name_ = other.problem_name_;
  problem_id_ = other.problem_id_;
  problem_type_ = other.problem_type_;

  evaluations_ = other.evaluations_;
  evaluations_constraints_ = other.evaluations_constraints_;
  final_target_delta_ = other.final_target_delta_;
  best_observed_fvalue_ = other.best_observed_fvalue_;
  best_observed_evaluation_ = other.best_observed_evaluation_;

  suite_ = other.suite_;
  suite_dep_index_ = other.suite_dep_index_;
  suite_dep_function_ = other.suite_dep_function_;
  suite_dep_instance_ = other.suite_dep_instance_;
}

void Problem::CountEvaluation(std::span<const double> x, std::span<const double> y) {
  evaluations_++; /* each derived class has its own counter, only the most outer will be visible */

  /* A little bit of bookkeeping */
  if (y[0] < best_observed_fvalue_) {
    bool feasible = true;
    if (number_of_constraints_ > 0) {
      std::vector<double> z(number_of_constraints_);
      feasible = IsFeasible(x, z);
    }
    if (feasible) {
      best_observed_fvalue_ = y[0];
      best_observed_evaluation_ = evaluations_;
    }
  }
}

void Problem::SetSuiteDependent(Suite *suite, const size_t index, const size_t function, const size_t instance) {
  this->suite_ = suite;
  suite_dep_index_ = index;
  suite_dep_function_ = function;
  suite_dep_instance_ = instance;
}

/**
 * @brief Initializes a single-objective problem using scalar values for smallest_value_of_interest,
 * largest_value_of_interest and best_parameter. Assumes all variables are continuous.
 */
void Problem::InitFromScalars(std::string_view name, const double smallest, const double largest, const double best) {

  assert(number_of_objectives_ == 1);
  assert(number_of_constraints_ == 0);
  problem_name_ = name;

  for (size_t i = 0; i < number_of_variables_; ++i) {
    smallest_values_of_interest_[i] = smallest;
    largest_values_of_interest_[i] = largest;
    best_parameter_[i] = best;
  }
  number_of_integer_variables_ = 0;
}

/***********************************************************************************************************/

bool Problem::IdIsFine(const std::string &id) {
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
void Problem::SetId(const char *format, ...) {
  va_list args;

  va_start(args, format);
  problem_id_ = StringVFormat(format, args);
  va_end(args);
  if (!IdIsFine(problem_id_)) {
    coco_error("Problem id should only contain standard chars, not like '%s'", problem_id_.c_str());
  }
}

/**
 * @brief Sets the problem_name using formatted printing (as in printf).
 */
void Problem::SetName(const char *format, ...) {
  va_list args;

  va_start(args, format);
  problem_name_ = StringVFormat(format, args);
  va_end(args);
}

/**
 * @brief Sets the problem_type using formatted printing (as in printf).
 */
void Problem::SetType(const char *format, ...) {
  va_list args;

  va_start(args, format);
  problem_type_ = StringVFormat(format, args);
  va_end(args);
}

/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods regarding the transformed COCO problem
 */
/**@{*/

/**
 * @brief Constructs a transformed problem that wraps the inner_problem.
 *
 * By default all methods will dispatch to the inner_problem. A prefix is prepended to the problem name
 * in order to reflect the transformation somewhere.
 */
TransformedProblem::TransformedProblem(ProblemPtr inner, std::string_view name_prefix)
    : Problem(inner->number_of_variables(), inner->number_of_objectives(), inner->number_of_constraints()),
      inner_problem_(std::move(inner)) {
  CopyPropertiesFrom(*inner_problem_);
  SetName("%.*s(%s)", static_cast<int>(name_prefix.size()), name_prefix.data(), inner_problem_->name().c_str());
}

/**
 * @brief Evaluates the objective function of the inner problem.
 */
void TransformedProblem::EvaluateFunction(std::span<const double> x, std::span<double> y) {
  inner_problem_->EvaluateFunction(x, y);
}

/**
 * @brief Evaluates the constraints of the inner problem.
 */
void TransformedProblem::EvaluateConstraint(std::span<const double> x, std::span<double> y, int update_counter) {
  inner_problem_->EvaluateConstraint(x, y, update_counter);
}

/**
 * @brief Calls the bbob_evaluate_gradient function on the inner problem.
 */
void TransformedProblem::EvaluateGradient(std::span<const double> x, std::span<double> y) {
  inner_problem_->EvaluateGradient(x, y);
}

/**
 * @brief Calls the coco_recommend_solution function on the inner problem.
 */
void TransformedProblem::RecommendSolution(std::span<const double> x) { inner_problem_->RecommendSolution(x); }

/**
 * @brief Maps base_x to the coordinates of the inner problem and from there to the coordinates of this problem.
 */
bool TransformedProblem::FromBaseCoordinates(std::span<const double> base_x, std::span<double> x) const {
  std::vector<double> inner_x(inner_problem_->number_of_variables());
  return inner_problem_->FromBaseCoordinates(base_x, inner_x) && Inverse(inner_x, x);
}

void TransformedProblem::UpdateBestParameter() {
  if (best_parameter_.empty())
    return;
  if (!Inverse(inner_problem_->best_parameter(), best_parameter_)) {
    coco_warning("%s: the transformation is not invertible, 'best_parameter' set to NAN", name().c_str());
    VectorSetToNan(best_parameter_);
  }
}

bool TransformedProblem::IdentityInverse(std::span<const double> inner_x, std::span<double> x) const {
  assert(inner_x.size() == x.size());
  std::copy(inner_x.begin(), inner_x.end(), x.begin());
  return true;
}

/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods regarding the bi-objective COCO problem
 */
/**@{*/

BiobjProblem::BiobjProblem(ProblemPtr objective1, ProblemPtr objective2,
                           std::span<const double> smallest_values_of_interest,
                           std::span<const double> largest_values_of_interest)
    : Problem(objective1->number_of_variables(), 2, 0), objective1_(std::move(objective1)),
      objective2_(std::move(objective2)) {
  assert(objective1_->number_of_objectives() == 1 && objective1_->number_of_constraints() == 0);
  assert(objective2_->number_of_objectives() == 1 && objective2_->number_of_constraints() == 0);
  assert(objective1_->number_of_variables() == objective2_->number_of_variables());

  SetId("%s__%s", objective1_->id().c_str(), objective2_->id().c_str());
  SetName("%s", id().c_str());
  SetType("%s_%s", objective1_->type().c_str(), objective2_->type().c_str());

  for (size_t i = 0; i < number_of_variables_; ++i)
    SetRegionOfInterest(i, smallest_values_of_interest[i], largest_values_of_interest[i]);
  assert(objective1_->number_of_integer_variables() == objective2_->number_of_integer_variables());
  number_of_integer_variables_ = objective1_->number_of_integer_variables();

  /* Compute the ideal and nadir values */
  best_value_[0] = objective1_->best_value()[0];
  best_value_[1] = objective2_->best_value()[0];
  coco_evaluate_function(objective1_.get(), objective2_->best_parameter().data(), &nadir_value_[0]);
  coco_evaluate_function(objective2_.get(), objective1_->best_parameter().data(), &nadir_value_[1]);
}

/**
 * @brief Evaluates both objective functions.
 */
void BiobjProblem::EvaluateFunction(std::span<const double> x, std::span<double> y) {
  objective1_->EvaluateFunction(x, y.first(1));
  objective2_->EvaluateFunction(x, y.subspan(1));
}

/**@}*/

/***********************************************************************************************************/

} // namespace coco

/* C interface (declared in coco.h) */

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
  assert(problem != nullptr);
  FromC(problem)->EvaluateFunctionChecked({x, FromC(problem)->number_of_variables()},
                                          {y, FromC(problem)->number_of_objectives()});
}

void coco_evaluate_constraint(coco_problem_t *problem, const double *x, double *y) {
  assert(problem != nullptr);
  FromC(problem)->EvaluateConstraintChecked({x, FromC(problem)->number_of_variables()},
                                            {y, FromC(problem)->number_of_constraints()}, 1);
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
  FromC(problem)->RecommendSolution({x, coco_problem_get_dimension(problem)});
}

void coco_problem_free(coco_problem_t *problem) {
  assert(problem != nullptr);
  delete FromC(problem);
}

size_t coco_problem_get_evaluations(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return FromC(problem)->evaluations();
}

size_t coco_problem_get_evaluations_constraints(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return FromC(problem)->evaluations_constraints();
}

/**
 * @note Can be used to prevent unnecessary burning of CPU time.
 */
int coco_problem_final_target_hit(const coco_problem_t *problem) {
  assert(problem != nullptr);
  if (coco_problem_get_number_of_objectives(problem) != 1 || coco_problem_get_evaluations(problem) < 1)
    return 0;
  if (FromC(problem)->best_value().empty())
    return 0;
  return FromC(problem)->best_observed_fvalue() <=
                 FromC(problem)->best_value()[0] + FromC(problem)->final_target_delta()
             ? 1
             : 0;
}

/**
 * @note Tentative...
 */
double coco_problem_get_best_observed_fvalue1(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return FromC(problem)->best_observed_fvalue();
}

/**
 * @brief Returns the optimal function value of the problem
 */
double coco_problem_get_best_value(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(!FromC(problem)->best_value().empty());
  return FromC(problem)->best_value()[0];
}

/**
 * @note This function breaks the black-box property: the returned  value is not
 * meant to be used by the optimization algorithm.
 */
double coco_problem_get_final_target_fvalue1(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(!FromC(problem)->best_value().empty());
  return FromC(problem)->best_value()[0] + FromC(problem)->final_target_delta();
}

/**
 * @note Do not modify the returned string! If you free the problem, the returned pointer becomes invalid.
 * When in doubt, copy the returned value.
 */
const char *coco_problem_get_name(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(!FromC(problem)->name().empty());
  return FromC(problem)->name().c_str();
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
  assert(!FromC(problem)->id().empty());
  return FromC(problem)->id().c_str();
}

const char *coco_problem_get_type(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(!FromC(problem)->type().empty());
  return FromC(problem)->type().c_str();
}

size_t coco_problem_get_dimension(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(FromC(problem)->number_of_variables() > 0);
  return FromC(problem)->number_of_variables();
}

size_t coco_problem_get_number_of_objectives(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return FromC(problem)->number_of_objectives();
}

size_t coco_problem_get_number_of_constraints(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return FromC(problem)->number_of_constraints();
}

const double *coco_problem_get_smallest_values_of_interest(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(!FromC(problem)->smallest_values_of_interest().empty());
  return FromC(problem)->smallest_values_of_interest().data();
}

const double *coco_problem_get_largest_values_of_interest(const coco_problem_t *problem) {
  assert(problem != nullptr);
  assert(!FromC(problem)->largest_values_of_interest().empty());
  return FromC(problem)->largest_values_of_interest().data();
}

size_t coco_problem_get_number_of_integer_variables(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return FromC(problem)->number_of_integer_variables();
}

const double *coco_problem_get_largest_fvalues_of_interest(const coco_problem_t *problem) {
  assert(problem != nullptr);
  if (FromC(problem)->number_of_objectives() == 1)
    coco_error(
        "coco_problem_get_largest_fvalues_of_interest(): f-values of interest undefined for single-objective problems");
  if (FromC(problem)->nadir_value().empty())
    coco_error("coco_problem_get_largest_fvalues_of_interest(): f-values of interest undefined");
  return FromC(problem)->nadir_value().data();
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
  if (!FromC(problem)->initial_solution().empty()) {
    for (size_t i = 0; i < FromC(problem)->number_of_variables(); ++i)
      initial_solution[i] = FromC(problem)->initial_solution()[i];
  } else {
    assert(!FromC(problem)->smallest_values_of_interest().empty());
    assert(!FromC(problem)->largest_values_of_interest().empty());
    for (size_t i = 0; i < FromC(problem)->number_of_variables(); ++i)
      initial_solution[i] =
          FromC(problem)->smallest_values_of_interest()[i] +
          0.5 * (FromC(problem)->largest_values_of_interest()[i] - FromC(problem)->smallest_values_of_interest()[i]);
    if (FromC(problem)->number_of_integer_variables() > 0) {
      for (size_t i = 0; i < FromC(problem)->number_of_integer_variables(); ++i) {
        initial_solution[i] = coco::DoubleRound(initial_solution[i]);
      }
    }
  }
}

size_t coco_problem_get_suite_dep_index(const coco_problem_t *problem) {
  assert(problem != nullptr);
  return FromC(problem)->suite_dep_index();
}

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

  if (FromC(problem)->best_parameter().empty() || FromC(problem)->number_of_objectives() != 1) {
    return 0;
  }

  // No more benchmarking on this problem!
  FromC(problem)->Taint();

  for (size_t i = 0; i < FromC(problem)->number_of_variables(); ++i) {
    best_parameter[i] = FromC(problem)->best_parameter()[i];
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
int coco_problem_is_tainted(const coco_problem_t *problem) { return FromC(problem)->is_tainted(); }

void bbob_problem_best_parameter_print(const coco_problem_t *problem) {
  coco::File file;
  assert(problem != nullptr);
  assert(!FromC(problem)->best_parameter().empty());
  if (file.Open("._bbob_problem_best_parameter.txt", "w")) {
    for (size_t i = 0; i < FromC(problem)->number_of_variables(); ++i)
      fprintf(file, " %.16f ", FromC(problem)->best_parameter()[i]);
  }
}

void bbob_biobj_problem_best_parameter_print(const coco_problem_t *problem) {
  coco::File file;

  /* We need the underlying single objective problems of the bi-objective problem, which are not accessible
   * once the problem has been wrapped, e.g. with a logger.
   */
  const auto *biobj = dynamic_cast<const coco::BiobjProblem *>(FromC(problem));
  if (biobj == nullptr) {
    coco_warning("I'm sorry, Dave. I'm afraid I can't do that.");
    /* In case this function was called previously, overwrite any results. */
    file.Open("._bbob_biobj_problem_best_parameter.txt", "w");
    return;
  }

  const coco::Problem *problem1 = biobj->objective1();
  const coco::Problem *problem2 = biobj->objective2();
  assert(!problem1->best_parameter().empty());
  assert(!problem2->best_parameter().empty());

  if (file.Open("._bbob_biobj_problem_best_parameter.txt", "w")) {
    for (size_t i = 0; i < biobj->number_of_variables(); ++i)
      fprintf(file, " %.16f ", problem1->best_parameter()[i]);
    fprintf(file, "\n");
    for (size_t i = 0; i < biobj->number_of_variables(); ++i)
      fprintf(file, " %.16f ", problem2->best_parameter()[i]);
  }
}
