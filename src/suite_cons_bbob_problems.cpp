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
 *       objective function, e.g. EllipsoidConsBbobProblemAllocate().
 *
 * @note Steps 4 and 5 are done within CLinearConsBbobProblemAllocate().
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
 * First, we create a Problem object for the objective function.
 * Then, Problem objects are created for each constraint and
 * stacked together into a single Problem object (see "c_linear.cpp").
 * Finally, the Problem object containing the constraints and
 * the Problem object containing the objective function are
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

namespace coco {

/**
 * @brief Maps the initial solution of the base problem through the inverse of the transformations on top of it.
 *
 * If the result plus xopt is not within (-5, 5) in every coordinate, the initial solution of the base problem
 * is halved until it is (xopt is needed because transform_vars_shift is applied afterwards). The initial
 * solution is left unchanged if the result contains NaN.
 */
static void ConsBbobInverseInitialSolution(Problem *problem, std::span<const double> xopt) {
  const size_t n = problem->number_of_variables();
  const double halving_factor = .5;
  std::vector<double> halved(n), inverted(n), sol(n);
  bool is_in_bounds = true;
  double di = 0.0;

  for (size_t j = 0;; ++j) {
    for (size_t i = 0; i < n; ++i)
      halved[i] = problem->initial_solution()[i] * pow(halving_factor, static_cast<double>(j));
    if (!problem->FromBaseCoordinates(halved, inverted))
      coco_error("ConsBbobInverseInitialSolution(): the transformations of %s are not invertible",
                 problem->id().c_str());
    for (size_t i = 0; i < n; ++i) {
      di = inverted[i];
      const double xi = di + xopt[i];
      is_in_bounds = (-5.0 < xi && xi < 5.0);
      if (!is_in_bounds)
        break;
      sol[i] = di;
    }
    if (is_in_bounds || IsNan(di))
      break;
  }
  if (!VectorContainsNan(sol))
    problem->SetInitialSolution(sol);
}

/**
 * @brief Calculates the objective function type based on the value
 *        of "function"
 */
static size_t ObjFunctionType(const size_t function) {

  size_t problems_per_obj_function_type = 6;
  return static_cast<size_t>(ceil(static_cast<double>(function) / static_cast<double>(problems_per_obj_function_type)));
}

/**
 * @brief Returns the number of linear constraints associated to the
 *        value of "function"
 */
static size_t NbOfLinearConstraints(const size_t function, const size_t dimension) {

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
static void FeasibleDirectionSetLength(std::span<double> feasible_direction, std::span<const double> xopt, long rseed) {
  const size_t dimension = feasible_direction.size();
  assert(xopt.size() == dimension);
  const long seed_offset = 412;    /* was sampled uniform in 0-999 */
  const double feas_shrink = 0.75; /* scale randomly between 0.75 and 1.0 */
  const double feas_bound = 5.0;

  double r[1], maxabs, maxrel;

  maxabs = maxrel = 0.0;
  for (size_t i = 0; i < dimension; ++i) {
    maxabs = DoubleMax(maxabs, fabs(xopt[i]));
    maxrel = DoubleMax(maxrel, feasible_direction[i] / (feas_bound - xopt[i]));
    maxrel = DoubleMax(maxrel, feasible_direction[i] / (-feas_bound - xopt[i]));
  }
  if (maxabs > 4.01)
    coco_warning("feasible_direction_set_length: a component of fabs(xopt) was greater than 4.01");
  if (maxabs > 5.0)
    coco_error("feasible_direction_set_length: a component of fabs(xopt) was greater than 5.0");
  Bbob2009Unif(r, rseed + seed_offset);
  VectorScale(feasible_direction, feas_shrink + r[0] * (1 - feas_shrink), /* nominator */
              maxrel);                                                    /* denominator */
}

/**
 * @brief Objective function: sphere
 *        Constraint(s): linear
 */
static ProblemPtr SphereCLinearConsBbobProblemAllocate(const size_t function, const size_t dimension,
                                                       const size_t instance, const size_t number_of_linear_constraints,
                                                       const long rseed, std::span<double> feasible_direction,
                                                       std::span<const double> xopt, const char *problem_id_template,
                                                       const char *problem_name_template) {

  ProblemPtr problem;
  ProblemPtr problem_c;

  const std::vector<double> all_zeros(dimension, 0.0);

  /* Create the objective function */
  problem =
      SphereConsBbobProblemAllocate(function, dimension, instance, rseed, problem_id_template, problem_name_template);

  problem->EvaluateGradient(all_zeros, feasible_direction);
  FeasibleDirectionSetLength(feasible_direction, xopt, rseed);

  /* Create the constraints. Use the gradient of the objective
   * function at the origin to build the first constraint.
   */
  problem_c = CLinearConsBbobProblemAllocate(function, dimension, instance, number_of_linear_constraints,
                                             problem_id_template, problem_name_template, feasible_direction);

  const std::string problem_type_temp = problem->type();
  const std::string problem_type_constraints = problem_c->type();
  problem = std::make_unique<StackedProblem>(std::move(problem), std::move(problem_c));

  /* The best parameter of the stacked problem is the origin, store its
   * objective function value into problem->best_value.
   */
  problem->EvaluateBestValue();

  /* To be sure that everything is correct, reset the f-evaluations
   * and g-evaluations counters to zero.
   */
  problem->ResetEvaluationCounters();

  /* Apply a translation to the whole problem so that the constrained
   * minimum is no longer at the origin.
   */
  problem = std::make_unique<TransformVarsShift>(std::move(problem), xopt, 0);

  /* Construct problem type */
  problem->SetType("%s_%s", problem_type_temp.c_str(), problem_type_constraints.c_str());

  return problem;
}

/**
 * @brief Objective function: ellipsoid
 *        Constraint(s): linear
 */
static ProblemPtr EllipsoidCLinearConsBbobProblemAllocate(const size_t function, const size_t dimension,
                                                          const size_t instance,
                                                          const size_t number_of_linear_constraints, const long rseed,
                                                          std::span<double> feasible_direction,
                                                          std::span<const double> xopt, const char *problem_id_template,
                                                          const char *problem_name_template) {

  ProblemPtr problem;
  ProblemPtr problem_c;

  const std::vector<double> all_zeros(dimension, 0.0);

  /* Create the objective function */
  problem = EllipsoidConsBbobProblemAllocate(function, dimension, instance, rseed, BbobVariant::BBOB,
                                             problem_id_template, problem_name_template);

  problem->EvaluateGradient(all_zeros, feasible_direction);
  FeasibleDirectionSetLength(feasible_direction, xopt, rseed);

  /* Create the constraints. Use the gradient of the objective
   * function at the origin to build the first constraint.
   */
  problem_c = CLinearConsBbobProblemAllocate(function, dimension, instance, number_of_linear_constraints,
                                             problem_id_template, problem_name_template, feasible_direction);

  const std::string problem_type_temp = problem->type();

  /* Build the final constrained function by stacking the objective
   * function and the constraints into one Problem type.
   */
  const std::string problem_type_constraints = problem_c->type();
  problem = std::make_unique<StackedProblem>(std::move(problem), std::move(problem_c));

  /* The best parameter of the stacked problem is the origin, store its
   * objective function value into problem->best_value.
   */
  problem->EvaluateBestValue();

  /* To be sure that everything is correct, reset the f-evaluations
   * and g-evaluations counters to zero.
   */
  problem->ResetEvaluationCounters();

  problem = std::make_unique<TransformVarsOscillate>(std::move(problem));
  ConsBbobInverseInitialSolution(problem.get(), xopt);

  /* Apply a translation to the whole problem so that the constrained
   * minimum is no longer at the origin.
   */
  problem = std::make_unique<TransformVarsShift>(std::move(problem), xopt, 0);

  /* Construct problem type */
  problem->SetType("%s_%s", problem_type_temp.c_str(), problem_type_constraints.c_str());

  return problem;
}

/**
 * @brief Objective function: rotated ellipsoid
 *        Constraint(s): linear
 */
static ProblemPtr
EllipsoidRotatedCLinearConsBbobProblemAllocate(const size_t function, const size_t dimension, const size_t instance,
                                               const size_t number_of_linear_constraints, const long rseed,
                                               std::span<double> feasible_direction, std::span<const double> xopt,
                                               const char *problem_id_template, const char *problem_name_template) {

  ProblemPtr problem;
  ProblemPtr problem_c;

  const std::vector<double> all_zeros(dimension, 0.0);

  /* Create the objective function */
  problem = EllipsoidRotatedConsBbobProblemAllocate(function, dimension, instance, rseed, BbobVariant::BBOB,
                                                    problem_id_template, problem_name_template);

  problem->EvaluateGradient(all_zeros, feasible_direction);
  FeasibleDirectionSetLength(feasible_direction, xopt, rseed);

  /* Create the constraints. Use the gradient of the objective
   * function at the origin to build the first constraint.
   */
  problem_c = CLinearConsBbobProblemAllocate(function, dimension, instance, number_of_linear_constraints,
                                             problem_id_template, problem_name_template, feasible_direction);

  const std::string problem_type_temp = problem->type();

  /* Build the final constrained function by stacking the objective
   * function and the constraints into one Problem type.
   */
  const std::string problem_type_constraints = problem_c->type();
  problem = std::make_unique<StackedProblem>(std::move(problem), std::move(problem_c));

  /* The best parameter of the stacked problem is the origin, store its
   * objective function value into problem->best_value.
   */
  problem->EvaluateBestValue();

  /* To be sure that everything is correct, reset the f-evaluations
   * and g-evaluations counters to zero.
   */
  problem->ResetEvaluationCounters();

  problem = std::make_unique<TransformVarsOscillate>(std::move(problem));
  ConsBbobInverseInitialSolution(problem.get(), xopt);

  /* Apply a translation to the whole problem so that the constrained
   * minimum is no longer at the origin .
   */
  problem = std::make_unique<TransformVarsShift>(std::move(problem), xopt, 0);

  /* Construct problem type */
  problem->SetType("%s_%s", problem_type_temp.c_str(), problem_type_constraints.c_str());

  return problem;
}

/**
 * @brief Objective function: linear slope
 *        Constraint(s): linear
 */
static ProblemPtr
LinearSlopeCLinearConsBbobProblemAllocate(const size_t function, const size_t dimension, const size_t instance,
                                          const size_t number_of_linear_constraints, const long rseed,
                                          std::span<double> feasible_direction, std::span<const double> xopt,
                                          const char *problem_id_template, const char *problem_name_template) {

  ProblemPtr problem;
  ProblemPtr problem_c;

  const std::vector<double> all_zeros(dimension, 0.0);

  /* Create the objective function */
  problem = LinearSlopeConsBbobProblemAllocate(function, dimension, instance, rseed, problem_id_template,
                                               problem_name_template);

  problem->EvaluateGradient(all_zeros, feasible_direction);
  FeasibleDirectionSetLength(feasible_direction, xopt, rseed);

  /* Create the constraints. Use the gradient of the objective
   * function at the origin to build the first constraint.
   */
  problem_c = CLinearConsBbobProblemAllocate(function, dimension, instance, number_of_linear_constraints,
                                             problem_id_template, problem_name_template, feasible_direction);

  const std::string problem_type_temp = problem->type();

  /* Build the final constrained function by stacking the objective
   * function and the constraints into one Problem type.
   */
  const std::string problem_type_constraints = problem_c->type();
  problem = std::make_unique<StackedProblem>(std::move(problem), std::move(problem_c));

  /* The best parameter of the stacked problem is the origin, store its
   * objective function value into problem->best_value.
   */

  /* Apply a translation to the whole problem so that the constrained
   * minimum is no longer at the origin.
   */
  problem = std::make_unique<TransformVarsShift>(std::move(problem), xopt, 1);

  problem->EvaluateBestValue();

  /* To be sure that everything is correct, reset the f-evaluations
   * and g-evaluations counters to zero.
   */
  problem->ResetEvaluationCounters();

  /* Construct problem type */
  problem->SetType("%s_%s", problem_type_temp.c_str(), problem_type_constraints.c_str());

  return problem;
}

/**
 * @brief Objective function: discus
 *        Constraint(s): linear
 */
static ProblemPtr DiscusCLinearConsBbobProblemAllocate(const size_t function, const size_t dimension,
                                                       const size_t instance, const size_t number_of_linear_constraints,
                                                       const long rseed, std::span<double> feasible_direction,
                                                       std::span<const double> xopt, const char *problem_id_template,
                                                       const char *problem_name_template) {

  ProblemPtr problem;
  ProblemPtr problem_c;

  const std::vector<double> all_zeros(dimension, 0.0);

  /* Create the objective function */
  problem = DiscusConsBbobProblemAllocate(function, dimension, instance, rseed, BbobVariant::BBOB, problem_id_template,
                                          problem_name_template);

  problem->EvaluateGradient(all_zeros, feasible_direction);
  FeasibleDirectionSetLength(feasible_direction, xopt, rseed);

  /* Create the constraints. Use the gradient of the objective
   * function at the origin to build the first constraint.
   */
  problem_c = CLinearConsBbobProblemAllocate(function, dimension, instance, number_of_linear_constraints,
                                             problem_id_template, problem_name_template, feasible_direction);

  const std::string problem_type_temp = problem->type();

  /* Build the final constrained function by stacking the objective
   * function and the constraints into one Problem type.
   */
  const std::string problem_type_constraints = problem_c->type();
  problem = std::make_unique<StackedProblem>(std::move(problem), std::move(problem_c));

  /* The best parameter of the stacked problem is the origin, store its
   * objective function value into problem->best_value.
   */
  problem->EvaluateBestValue();

  /* To be sure that everything is correct, reset the f-evaluations
   * and g-evaluations counters to zero.
   */
  problem->ResetEvaluationCounters();

  problem = std::make_unique<TransformVarsOscillate>(std::move(problem));
  ConsBbobInverseInitialSolution(problem.get(), xopt);

  /* Apply a translation to the whole problem so that the constrained
   * minimum is no longer at the origin.
   */
  problem = std::make_unique<TransformVarsShift>(std::move(problem), xopt, 0);

  /* Construct problem type */
  problem->SetType("%s_%s", problem_type_temp.c_str(), problem_type_constraints.c_str());

  return problem;
}

/**
 * @brief Objective function: bent cigar
 *        Constraint(s): linear
 */
static ProblemPtr BentCigarCLinearConsBbobProblemAllocate(const size_t function, const size_t dimension,
                                                          const size_t instance,
                                                          const size_t number_of_linear_constraints, const long rseed,
                                                          std::span<double> feasible_direction,
                                                          std::span<const double> xopt, const char *problem_id_template,
                                                          const char *problem_name_template) {

  ProblemPtr problem;
  ProblemPtr problem_c;

  const std::vector<double> all_zeros(dimension, 0.0);

  /* Create the objective function */
  problem = BentCigarConsBbobProblemAllocate(function, dimension, instance, rseed, BbobVariant::BBOB,
                                             problem_id_template, problem_name_template);

  problem->EvaluateGradient(all_zeros, feasible_direction);
  FeasibleDirectionSetLength(feasible_direction, xopt, rseed);

  /* Create the constraints. Use the gradient of the objective
   * function at the origin to build the first constraint.
   */
  problem_c = CLinearConsBbobProblemAllocate(function, dimension, instance, number_of_linear_constraints,
                                             problem_id_template, problem_name_template, feasible_direction);

  const std::string problem_type_temp = problem->type();

  /* Build the final constrained function by stacking the objective
   * function and the constraints into one Problem type.
   */
  const std::string problem_type_constraints = problem_c->type();
  problem = std::make_unique<StackedProblem>(std::move(problem), std::move(problem_c));

  /* The best parameter of the stacked problem is the origin, store its
   * objective function value into problem->best_value.
   */
  problem->EvaluateBestValue();

  /* To be sure that everything is correct, reset the f-evaluations
   * and g-evaluations counters to zero.
   */
  problem->ResetEvaluationCounters();

  problem = std::make_unique<TransformVarsAsymmetric>(std::move(problem), 0.2);
  ConsBbobInverseInitialSolution(problem.get(), xopt);

  /* Apply a translation to the whole problem so that the constrained
   * minimum is no longer at the origin.
   */
  problem = std::make_unique<TransformVarsShift>(std::move(problem), xopt, 0);

  /* Construct problem type */
  problem->SetType("%s_%s", problem_type_temp.c_str(), problem_type_constraints.c_str());

  return problem;
}

/**
 * @brief Objective function: different powers
 *        Constraint(s): linear
 */
static ProblemPtr
DifferentPowersCLinearConsBbobProblemAllocate(const size_t function, const size_t dimension, const size_t instance,
                                              const size_t number_of_linear_constraints, const long rseed,
                                              std::span<double> feasible_direction, std::span<const double> xopt,
                                              const char *problem_id_template, const char *problem_name_template) {

  ProblemPtr problem;
  ProblemPtr problem_c;

  const std::vector<double> all_zeros(dimension, 0.0);

  /* Create the objective function */
  problem = DifferentPowersBbobConstrainedProblemAllocate(function, dimension, instance, rseed, BbobVariant::BBOB,
                                                          problem_id_template, problem_name_template);

  problem->EvaluateGradient(all_zeros, feasible_direction);
  FeasibleDirectionSetLength(feasible_direction, xopt, rseed);

  /* Create the constraints. Use the gradient of the objective
   * function at the origin to build the first constraint.
   */
  problem_c = CLinearConsBbobProblemAllocate(function, dimension, instance, number_of_linear_constraints,
                                             problem_id_template, problem_name_template, feasible_direction);

  const std::string problem_type_temp = problem->type();

  /* Build the final constrained function by stacking the objective
   * function and the constraints into one Problem type.
   */
  const std::string problem_type_constraints = problem_c->type();
  problem = std::make_unique<StackedProblem>(std::move(problem), std::move(problem_c));

  /* The best parameter of the stacked problem is the origin, store its
   * objective function value into problem->best_value.
   */
  problem->EvaluateBestValue();

  /* To be sure that everything is correct, reset the f-evaluations
   * and g-evaluations counters to zero.
   */
  problem->ResetEvaluationCounters();

  /* Apply a translation to the whole problem so that the constrained
   * minimum is no longer at the origin.
   */
  problem = std::make_unique<TransformVarsShift>(std::move(problem), xopt, 0);

  /* Construct problem type */
  problem->SetType("%s_%s", problem_type_temp.c_str(), problem_type_constraints.c_str());

  return problem;
}

/**
 * @brief Objective function: rastrigin
 *        Constraint(s): linear
 */
static ProblemPtr RastriginCLinearConsBbobProblemAllocate(const size_t function, const size_t dimension,
                                                          const size_t instance,
                                                          const size_t number_of_linear_constraints, const long rseed,
                                                          std::span<double> feasible_direction,
                                                          std::span<const double> xopt, const char *problem_id_template,
                                                          const char *problem_name_template) {

  ProblemPtr problem;
  ProblemPtr problem_c;

  const std::vector<double> all_zeros(dimension, 0.0);

  /* Create the objective function */
  problem = RastriginConsBbobProblemAllocate(function, dimension, instance, rseed, problem_id_template,
                                             problem_name_template);

  problem->EvaluateGradient(all_zeros, feasible_direction);
  FeasibleDirectionSetLength(feasible_direction, xopt, rseed);

  /* Create the constraints. Use the feasible direction above
   * to build the first constraint.
   */
  problem_c = CLinearConsBbobProblemAllocate(function, dimension, instance, number_of_linear_constraints,
                                             problem_id_template, problem_name_template, feasible_direction);

  const std::string problem_type_temp = problem->type();

  /* Build the final constrained function by stacking the objective
   * function and the constraints into one Problem type.
   */
  const std::string problem_type_constraints = problem_c->type();
  problem = std::make_unique<StackedProblem>(std::move(problem), std::move(problem_c));

  /* The best parameter of the stacked problem is the origin, store its
   * objective function value into problem->best_value.
   */
  problem->EvaluateBestValue();

  /* To be sure that everything is correct, reset the f-evaluations
   * and g-evaluations counters to zero.
   */
  problem->ResetEvaluationCounters();

  problem = std::make_unique<TransformVarsAsymmetric>(std::move(problem), 0.2);
  problem = std::make_unique<TransformVarsOscillate>(std::move(problem));
  ConsBbobInverseInitialSolution(problem.get(), xopt);

  /* Apply a translation to the whole problem so that the constrained
   * minimum is no longer at the origin.
   */
  problem = std::make_unique<TransformVarsShift>(std::move(problem), xopt, 0);

  /* Construct problem type */
  problem->SetType("%s_%s", problem_type_temp.c_str(), problem_type_constraints.c_str());

  return problem;
}

/**
 * @brief Objective function: rastrigin
 *        Constraint(s): linear
 */
static ProblemPtr
RastriginRotatedCLinearConsBbobProblemAllocate(const size_t function, const size_t dimension, const size_t instance,
                                               const size_t number_of_linear_constraints, const long rseed,
                                               std::span<double> feasible_direction, std::span<const double> xopt,
                                               const char *problem_id_template, const char *problem_name_template) {

  ProblemPtr problem;
  ProblemPtr problem_c;

  const std::vector<double> all_zeros(dimension, 0.0);

  /* Create the objective function */
  problem = RastriginRotatedConsBbobProblemAllocate(function, dimension, instance, rseed, problem_id_template,
                                                    problem_name_template);

  problem->EvaluateGradient(all_zeros, feasible_direction);
  FeasibleDirectionSetLength(feasible_direction, xopt, rseed);

  /* Create the constraints. Use the feasible direction above
   * to build the first constraint.
   */
  problem_c = CLinearConsBbobProblemAllocate(function, dimension, instance, number_of_linear_constraints,
                                             problem_id_template, problem_name_template, feasible_direction);

  const std::string problem_type_temp = problem->type();

  /* Build the final constrained function by stacking the objective
   * function and the constraints into one Problem type.
   */
  const std::string problem_type_constraints = problem_c->type();
  problem = std::make_unique<StackedProblem>(std::move(problem), std::move(problem_c));

  /* The best parameter of the stacked problem is the origin, store its
   * objective function value into problem->best_value.
   */
  problem->EvaluateBestValue();

  /* To be sure that everything is correct, reset the f-evaluations
   * and g-evaluations counters to zero.
   */
  problem->ResetEvaluationCounters();

  problem = std::make_unique<TransformVarsAsymmetric>(std::move(problem), 0.2);
  problem = std::make_unique<TransformVarsOscillate>(std::move(problem));
  ConsBbobInverseInitialSolution(problem.get(), xopt);

  /* Apply a translation to the whole problem so that the constrained
   * minimum is no longer at the origin.
   */
  problem = std::make_unique<TransformVarsShift>(std::move(problem), xopt, 0);

  /* Construct problem type */
  problem->SetType("%s_%s", problem_type_temp.c_str(), problem_type_constraints.c_str());

  return problem;
}

} // namespace coco
