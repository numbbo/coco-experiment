/**
 * @file  c_linear.cpp
 * @brief Implements the linear constraints for the suite of
 *        constrained problems.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "coco.h"
#include "coco_internal.h"
#include "coco_problem.cpp"
namespace coco {

static ProblemPtr CGuaranteeFeasiblePoint(ProblemPtr problem, std::span<const double> feasible_point);

static ProblemPtr CLinearTransform(ProblemPtr inner_problem, std::span<const double> gradient,
                                   const double x_shift_factor);

static ProblemPtr CLinearSingleConsBbobProblemAllocate(const size_t function, const size_t dimension,
                                                       const size_t instance, const size_t constraint_number,
                                                       const double factor1, const char *problem_id_template,
                                                       const char *problem_name_template, std::span<double> gradient,
                                                       double x_shift_factor,
                                                       std::span<const double> feasible_direction);

static ProblemPtr CLinearConsBbobProblemAllocate(const size_t function, const size_t dimension, const size_t instance,
                                                 const size_t number_of_active_constraints,
                                                 const char *problem_id_template, const char *problem_name_template,
                                                 std::span<const double> feasible_direction);

/**
 * @brief A linear constraint with all-ones gradient.
 */
class CSumVariables final : public Problem {
public:
  explicit CSumVariables(const size_t number_of_variables) : Problem(number_of_variables, 0, 1) {
    SetId("%s_d%02zu", "linearconstraint", number_of_variables);
    for (size_t i = 0; i < number_of_variables; ++i) {
      smallest_values_of_interest_[i] = -5.0;
      largest_values_of_interest_[i] = 5.0;
    }
  }

  /**
   * @brief Evaluates the linear constraint with all-ones gradient at
   *        the point 'x' and stores the result into 'y'.
   */
  void EvaluateConstraint(std::span<const double> x, std::span<double> y, int) override {
    size_t i;

    assert(number_of_constraints_ == 1);
    y[0] = 0.0;
    for (i = 0; i < number_of_variables_; ++i)
      /* prevent that the optimal solution is infeasible due to loss of precision when shifted */
      if (x[i] > 1e-11 || x[i] < 0)
        y[0] += x[i];
  }
};

/**
 * @brief A linear constraint whose gradient is given and whose feasible domain is extended by a shift factor.
 */
class CLinear final : public TransformedProblem {
public:
  CLinear(ProblemPtr inner, std::span<const double> gradient, const double x_shift_factor, const double gradient_norm)
      : TransformedProblem(std::move(inner), "gradient_linear_constraint"), gradient_(gradient.begin(), gradient.end()),
        transformed_x_(inner_problem_->number_of_variables()), x_shift_factor_(x_shift_factor),
        gradient_norm_(gradient_norm) {}

  /**
   * @brief Evaluates the linear constraint at the point 'x' and stores
   *        the result in 'y'.
   */
  void EvaluateConstraint(std::span<const double> x, std::span<double> y, int update_counter) override {
    size_t i;
    double factor;

    assert(number_of_constraints_ == 1);
    for (i = 0, factor = x_shift_factor_ / gradient_norm_; i < number_of_variables_; ++i)
      transformed_x_[i] = (gradient_[i]) * (x[i] - factor * gradient_[i]);

    inner_problem_->EvaluateConstraint(transformed_x_, y, update_counter);
  }

  /** @brief Inverts the transformation of the constraint; not possible if a gradient component is 0. */
  bool Inverse(std::span<const double> inner_x, std::span<double> x) const override {
    const double factor = x_shift_factor_ / gradient_norm_;
    for (size_t i = 0; i < number_of_variables_; ++i) {
      if (gradient_[i] == 0.0)
        return false;
      x[i] = inner_x[i] / gradient_[i] + factor * gradient_[i];
    }
    return true;
  }

  const std::vector<double> &gradient() const { return gradient_; }
  double x_shift_factor() const { return x_shift_factor_; }

  /** @brief Turns the constraint around by negating its gradient. */
  void Flip() {
    for (size_t i = 0; i < number_of_variables_; ++i)
      gradient_[i] *= -1.0;
  }

  /**
   * @brief Disguises the gradient of the objective function in this constraint: subtracts a random multiple of
   * linear_combination (a positive combination of the other active constraints' gradients) from the gradient
   * and restores its norm, which does not change the feasibility of the initial solution.
   */
  void DisguiseGradient(std::span<const double> linear_combination, const double random_uniform) {
    const size_t dimension = number_of_variables_;
    const double norm = gradient_norm_;
    double fac = VectorScalarProduct(linear_combination, gradient_);
    if (fac < 0)
      coco_error("scalar product between first constraint and linear combination = %f < 0 should be > 0", fac);
    fac = fac <= 0 ? static_cast<double>(dimension) : DoubleMin(static_cast<double>(dimension), norm * norm / fac);
    fac *= (1 + random_uniform) / (2 + 1. / static_cast<double>(dimension)); /* bound to <1 and >1/3 */
    for (size_t i = 0; i < dimension; ++i)
      gradient_[i] -= fac * linear_combination[i];
    fac = VectorNorm(gradient_);      /* new gradient norm */
    if (fac == 0)                     /* make sure that gradient did not become zero */
      coco_error("new vector norm = %f == 0, should be > 0", fac);
    VectorScale(gradient_, norm, fac); /* restore original norm */
  }

private:
  std::vector<double> gradient_;
  std::vector<double> transformed_x_;
  const double x_shift_factor_; /* shift solution by - factor * gradient */
  const double gradient_norm_;
};

/**
 * @brief Returns the linear constraint behind the given problem pointer.
 */
static CLinear *CLinearGet(const Problem *problem) {
  assert(dynamic_cast<const CLinear *>(problem) != nullptr);
  return static_cast<CLinear *>(const_cast<Problem *>(problem));
}

/**
 * @brief Guarantees that "feasible_direction" is feasible w.r.t.
 *        the constraint in "problem" and records it as the
 *        initial feasible solution to this coco_problem.
 */
static ProblemPtr CGuaranteeFeasiblePoint(ProblemPtr problem, std::span<const double> feasible_direction) {

  CLinear *data;
  double constraint_value = 0.0;

  data = CLinearGet(problem.get());

  assert(problem->number_of_constraints() == 1);

  /* Let p be the gradient of the constraint in "problem".
   * Check whether p' * (feasible_direction) <= 0.
   */
  problem->EvaluateConstraint(feasible_direction, {&constraint_value, 1}, 0);

  /* Flip the constraint in "problem" if feasible_direction
   * is not feasible w.r.t. the constraint in "problem".
   */
  if (constraint_value > 0)
    data->Flip();

  problem->SetInitialSolution(feasible_direction);

  data = nullptr;
  return problem;
}

/**
 * @brief Transforms a linear constraint with all-ones gradient
 *        into a linear constraint whose gradient is passed
 *        as argument and with the feasible domain extended by
 *        shift_factor.
 */
static ProblemPtr CLinearTransform(ProblemPtr inner_problem, std::span<const double> gradient,
                                   const double shift_factor) {

  assert(gradient.size() == inner_problem->number_of_variables());
  double gradient_norm = VectorNorm(gradient);

  if (gradient_norm <= 0)
    coco_error("CLinearTransform(): gradient norm %f<=0 zero", gradient_norm);

  return std::make_unique<CLinear>(std::move(inner_problem), gradient, shift_factor, gradient_norm);
}

/**
 * @brief Builds a Problem containing one single linear constraint.
 *
 * This function is called by CLinearConsBbobProblemAllocate(),
 * the central function that stacks all the constraints built by
 * CLinearSingleConsBbobProblemAllocate() into one single
 * Problem object.
 */
static ProblemPtr CLinearSingleConsBbobProblemAllocate(const size_t function, const size_t dimension,
                                                       const size_t instance, const size_t constraint_number,
                                                       const double factor1, const char *problem_id_template,
                                                       const char *problem_name_template, std::span<double> gradient,
                                                       const double x_shift_factor,
                                                       std::span<const double> feasible_direction) {

  ProblemPtr problem;
  double factor2;

  problem = std::make_unique<CSumVariables>(dimension);

  long seed_cons_i = static_cast<long>(function + 10000 * instance + 50000 * constraint_number);
  RandomState random_generator(static_cast<uint32_t>(seed_cons_i));

  /* The constraints gradients are scaled with random numbers
   * 10**U[0,1] and 10**U_i[0,2], where U[a, b] is uniform in [a,b]
   * and only U_i is drawn for each constraint individually.
   * The random number 10**U[0,1] is given by the variable 'factor1'
   * while the random number 10**U_i[0,2] is calculated below and
   * stored as 'factor2'. (The exception is when the number of
   * constraints is n+1, in which case 'factor2' defines a random
   * number 10**U_i[0,1])
   */

  factor2 = pow(100.0, random_generator.Uniform());

  /* Set the gradient of the linear constraint if it is given.
   * This should be the case of the construction of the first
   * linear constraint only.
   */
  if (!gradient.empty()) {

    VectorScale(gradient, factor1 * factor2, VectorNorm(gradient));
    problem = CLinearTransform(std::move(problem), gradient, x_shift_factor);

  } else { /* Randomly generate the gradient of the linear constraint */

    std::vector<double> gradient_linear_constraint(dimension);

    /* Generate a pseudorandom vector with distribution N_i(0, I)
     * and scale it with 'factor1' and 'factor2' (see comments above)
     */
    for (size_t i = 0; i < dimension; ++i)
      gradient_linear_constraint[i] =
          factor1 * random_generator.Normal() * factor2 / sqrt(static_cast<double>(dimension));

    problem = CLinearTransform(std::move(problem), gradient_linear_constraint, x_shift_factor);
  }

  /* Guarantee that the vector feasible_point is feasible w.r.t. to
   * this constraint and set it as the initial solution.
   * The initial solution will be copied later to the constrained function
   * Problem object once the objective function and the constraint(s)
   * are stacked together in a StackedProblem.
   */
  if (!feasible_direction.empty())
    problem = CGuaranteeFeasiblePoint(std::move(problem), feasible_direction);

  problem->SetId(problem_id_template, function, instance, dimension);
  problem->SetName(problem_name_template, function, instance, dimension);
  problem->SetType("linear");
  return problem;
}

/**
 * @brief helper function to successively compute a linear combination of
 *        constraint gradients, add current gradient if weight is != 0.
 *
 */
static void ConUpdateLinearCombination(std::span<double> linear_combination, const Problem *problem, double weight) {
  const CLinear *data;

  data = CLinearGet(problem);
  if (data->gradient().empty()) {
    if (weight != 0) {
      coco_error("ConUpdateLinearCombination(): gradient of constraint was zero");
    } else {
      coco_warning("ConUpdateLinearCombination(): gradient of constraint was zero");
    }
  } else {
    if (data->x_shift_factor() != 0)
      coco_warning("Inactive constraint passed to update_linear_combination, x_shift_factor=%f",
                   data->x_shift_factor());
    if (weight == 0)
      return; /* nothing to add */
    if (weight < 0)
      coco_warning("con_update_linear_combination: weight=%f < 0, should be > 0", weight);
    for (size_t i = 0; i < problem->number_of_variables(); ++i)
      linear_combination[i] += weight * data->gradient()[i];
  }
}

/**
 * @brief Builds a Problem containing all the linear constraints
 *        by stacking them all.
 *
 * The constraints' gradients are randomly generated with distribution
 * 10**U[0,1] * N_i(0, I) * 10**U_i[0,2], where U[a, b] is uniform
 * in [a,b] and only U_i is drawn for each constraint individually.
 * The exception is the first constraint, whose gradient is given by
 * 10**U[0,1] * (-feasible_direction) * 10**U_i[0,2].
 *
 * Each constraint is built by calling the function
 * CLinearSingleConsBbobProblemAllocate(), which returns a
 * Problem object that defines the constraint. The resulting
 * Problem objects are then stacked together into one single
 * Problem object that is returned by the function.
 */
static ProblemPtr CLinearConsBbobProblemAllocate(const size_t function, const size_t dimension, const size_t instance,
                                                 const size_t number_of_active_constraints,
                                                 const char *problem_id_template, const char *problem_name_template,
                                                 std::span<const double> feasible_direction) {

  const double global_scaling_factor = 100.;
  size_t i, j;

  ProblemPtr problem_c;
  ProblemPtr problem_c2;
  std::span<double> gradient;
  CLinear *first_constraint_data = nullptr;
  double shift_factor;
  double factor1;
  size_t inactive_constraints_left =
      number_of_active_constraints / 2; /* TODO: decide whether this should go into the interface */
  int disguise_gradient = 1;

  if (std::string_view(problem_id_template).substr(0, 28) == "bbob-constrained-active-only")
    inactive_constraints_left = 0;

  if (std::string_view(problem_id_template).substr(0, 28) == "bbob-constrained-no-disguise")
    disguise_gradient = 0;

  size_t number_of_linear_constraints = number_of_active_constraints + inactive_constraints_left;

  std::vector<double> gradient_c1(dimension);
  std::vector<double> linear_combination(dimension, 0.0);

  for (i = 0; i < dimension; ++i)
    gradient_c1[i] = -feasible_direction[i];

  /* Build a Problem object for each constraint.
   * The constraints' gradients are generated randomly with
   * distribution 10**U[0,1] * N_i(0, I/n) * 10**U_i[0,2]
   * where U[a, b] is uniform in [a,b] and only U_i is drawn
   * for each constraint individually.
   */

  /* Calculate the first random factor 10**U[0,1]. */
  long seed_cons = static_cast<long>(function + 10000 * instance);
  RandomState random_generator(static_cast<uint32_t>(seed_cons));
  RandomState random_generator2(static_cast<uint32_t>(seed_cons));
  factor1 = global_scaling_factor * pow(10.0, random_generator.Uniform());

  /* Build the first linear constraint using 'gradient_c1' to build
   * its gradient.
   */
  /* set gradient depending on instance number */
  gradient = instance % number_of_active_constraints ? std::span<double>() : std::span<double>(gradient_c1);
  problem_c = CLinearSingleConsBbobProblemAllocate(function, dimension, instance, 1, factor1, problem_id_template,
                                                   problem_name_template, gradient, 0.0, feasible_direction);
  if (!gradient.empty()) /* preserve the constraint based on function gradient */
    first_constraint_data = CLinearGet(problem_c.get());
  else
    ConUpdateLinearCombination(linear_combination, problem_c.get(), random_generator2.Uniform());

  /* Instantiate the other linear constraints (if any) and stack them
   * all into problem_c
   */
  for (j = 2, i = 1; j <= number_of_linear_constraints; ++j) {

    /* Instantiate a new problem containing one linear constraint only */
    /* set gradient depending on instance number */

    if (i < number_of_active_constraints &&
        random_generator.Uniform() + 1e-23 >
            static_cast<double>(inactive_constraints_left) /
                static_cast<double>(number_of_linear_constraints)) { /* create an active constraint */
      gradient =
          (++i - 1 + instance) % number_of_active_constraints ? std::span<double>() : std::span<double>(gradient_c1);
      shift_factor = 0.0;
    } else { /* create an inactive (shifted) constraint */
      if (inactive_constraints_left-- <= 0)
        coco_error("CLinearConsBbobProblemAllocate(): no inactive left i=%zu j=%zu nb_act=%zu nb_con=%zu", i, j,
                   number_of_active_constraints, number_of_linear_constraints);
      gradient = {};
      shift_factor = 0.01 + 2.0 * random_generator.Uniform();
    }

    problem_c2 = CLinearSingleConsBbobProblemAllocate(
        function, dimension, instance,
        i <= number_of_active_constraints ? i : number_of_linear_constraints - inactive_constraints_left, factor1,
        problem_id_template, problem_name_template, gradient, shift_factor, feasible_direction);
    if (shift_factor == 0) {  /* active constraint */
      if (!gradient.empty()) { /* preserve the constraint based on function gradient */
        if (first_constraint_data)
          coco_warning(
              "CLinearConsBbobProblemAllocate(): first_constraint_data already assigned, this is probably a bug");
        first_constraint_data = CLinearGet(problem_c2.get());
      } else
        ConUpdateLinearCombination(linear_combination, problem_c2.get(), random_generator2.Uniform());
    }

    const std::string problem_type_c2 = problem_c2->type();
    problem_c = std::make_unique<StackedProblem>(std::move(problem_c), std::move(problem_c2));

    /* Use the standard stacked problem_id as problem_name and
     * construct a new suite-specific problem_id
     */
    problem_c->SetName("%s", problem_c->id().c_str());
    problem_c->SetId("bbob-constrained_f%02zu_i%02zu_d%02zu", function, instance, dimension);

    /* Construct problem type */
    problem_c->SetType("%s_%s", problem_type_c2.c_str(), problem_type_c2.c_str());
  }

  /* Modify first constraint without changing the feasible solution,
   * thereby disguising the gradient of the function
   */
  if (number_of_active_constraints > 1 && disguise_gradient)
    first_constraint_data->DisguiseGradient(linear_combination, random_generator2.Uniform());

  return problem_c;
}

} // namespace coco
