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

/**
 * @brief A linear constraint whose feasible domain is extended by a shift factor.
 *
 * g(x) = sum_i t_i over all t_i = a_i * (x_i - s / |a| * a_i) with t_i > 1e-11 or t_i < 0, where a is the
 * gradient and s the shift factor. Leaving out the small positive t_i prevents that the optimal solution is
 * infeasible due to loss of precision when shifted.
 */
struct LinearConstraint {
  std::vector<double> gradient;
  double x_shift_factor; /* shift solution by - factor * gradient */
  double gradient_norm;

  double operator()(std::span<const double> x) const {
    const double factor = x_shift_factor / gradient_norm;
    double y = 0.0;
    for (size_t i = 0; i < gradient.size(); ++i) {
      const double t = gradient[i] * (x[i] - factor * gradient[i]);
      if (t > 1e-11 || t < 0)
        y += t;
    }
    return y;
  }

  /** @brief Turns the constraint around by negating its gradient. */
  void Flip() {
    for (double &gi : gradient)
      gi *= -1.0;
  }

  /**
   * @brief Disguises the gradient of the objective function in this constraint: subtracts a random multiple of
   * linear_combination (a positive combination of the other active constraints' gradients) from the gradient
   * and restores its norm, which does not change the feasibility of the initial solution.
   */
  void DisguiseGradient(std::span<const double> linear_combination, const double random_uniform) {
    const size_t dimension = gradient.size();
    const double norm = gradient_norm;
    double fac = VectorScalarProduct(linear_combination, gradient);
    if (fac < 0)
      coco_error("scalar product between first constraint and linear combination = %f < 0 should be > 0", fac);
    fac = fac <= 0 ? static_cast<double>(dimension) : DoubleMin(static_cast<double>(dimension), norm * norm / fac);
    fac *= (1 + random_uniform) / (2 + 1. / static_cast<double>(dimension)); /* bound to <1 and >1/3 */
    for (size_t i = 0; i < dimension; ++i)
      gradient[i] -= fac * linear_combination[i];
    fac = VectorNorm(gradient); /* new gradient norm */
    if (fac == 0)               /* make sure that gradient did not become zero */
      coco_error("new vector norm = %f == 0, should be > 0", fac);
    VectorScale(gradient, norm, fac); /* restore original norm */
  }
};

/**
 * @brief Builds one linear constraint.
 *
 * The gradient is given for the first active constraint only and drawn randomly otherwise. The constraint is
 * flipped if necessary so that feasible_direction is feasible.
 */
static LinearConstraint CLinearSingleConstraint(const size_t function, const size_t dimension, const size_t instance,
                                                const size_t constraint_number, const double factor1,
                                                std::span<double> gradient, const double x_shift_factor,
                                                std::span<const double> feasible_direction) {

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

  const double factor2 = pow(100.0, random_generator.Uniform());

  LinearConstraint constraint;

  /* Set the gradient of the linear constraint if it is given.
   * This should be the case of the construction of the first
   * linear constraint only.
   */
  if (!gradient.empty()) {

    VectorScale(gradient, factor1 * factor2, VectorNorm(gradient));
    constraint.gradient.assign(gradient.begin(), gradient.end());

  } else { /* Randomly generate the gradient of the linear constraint */

    constraint.gradient.resize(dimension);

    /* Generate a pseudorandom vector with distribution N_i(0, I)
     * and scale it with 'factor1' and 'factor2' (see comments above)
     */
    for (size_t i = 0; i < dimension; ++i)
      constraint.gradient[i] = factor1 * random_generator.Normal() * factor2 / sqrt(static_cast<double>(dimension));
  }

  constraint.x_shift_factor = x_shift_factor;
  constraint.gradient_norm = VectorNorm(constraint.gradient);
  if (constraint.gradient_norm <= 0)
    coco_error("CLinearSingleConstraint(): gradient norm %f<=0 zero", constraint.gradient_norm);

  /* Guarantee that feasible_direction is feasible w.r.t. this constraint */
  if (constraint(feasible_direction) > 0)
    constraint.Flip();

  return constraint;
}

/**
 * @brief helper function to successively compute a linear combination of
 *        constraint gradients, add current gradient if weight is != 0.
 *
 */
static void ConUpdateLinearCombination(std::span<double> linear_combination, const LinearConstraint &constraint,
                                       double weight) {
  if (constraint.x_shift_factor != 0)
    coco_warning("Inactive constraint passed to update_linear_combination, x_shift_factor=%f",
                 constraint.x_shift_factor);
  if (weight == 0)
    return; /* nothing to add */
  if (weight < 0)
    coco_warning("con_update_linear_combination: weight=%f < 0, should be > 0", weight);
  for (size_t i = 0; i < linear_combination.size(); ++i)
    linear_combination[i] += weight * constraint.gradient[i];
}

/**
 * @brief Builds the linear constraints of a constrained bbob problem.
 *
 * The constraints' gradients are randomly generated with distribution
 * 10**U[0,1] * N_i(0, I) * 10**U_i[0,2], where U[a, b] is uniform
 * in [a,b] and only U_i is drawn for each constraint individually.
 * The exception is the first constraint, whose gradient is given by
 * 10**U[0,1] * (-feasible_direction) * 10**U_i[0,2].
 *
 * Each constraint is built by calling the function CLinearSingleConstraint(). Besides the active constraints,
 * number_of_active_constraints / 2 inactive (shifted) ones are added unless active_only is set. Unless
 * disguise is unset, the gradient of the constraint based on the function gradient is disguised.
 */
static std::vector<LinearConstraint> CLinearConsBbobConstraints(const size_t function, const size_t dimension,
                                                                const size_t instance,
                                                                const size_t number_of_active_constraints,
                                                                const bool active_only, const bool disguise,
                                                                std::span<const double> feasible_direction) {

  const double global_scaling_factor = 100.;
  size_t i, j;

  std::vector<LinearConstraint> constraints;
  std::span<double> gradient;
  size_t first_constraint = SIZE_MAX; /* index of the constraint based on the function gradient */
  double shift_factor;
  double factor1;
  size_t inactive_constraints_left =
      active_only ? 0 : number_of_active_constraints / 2; /* TODO: decide whether this should go into the interface */

  size_t number_of_linear_constraints = number_of_active_constraints + inactive_constraints_left;

  std::vector<double> gradient_c1(dimension);
  std::vector<double> linear_combination(dimension, 0.0);

  for (i = 0; i < dimension; ++i)
    gradient_c1[i] = -feasible_direction[i];

  /* Build each constraint.
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
  constraints.push_back(
      CLinearSingleConstraint(function, dimension, instance, 1, factor1, gradient, 0.0, feasible_direction));
  if (!gradient.empty()) /* preserve the constraint based on function gradient */
    first_constraint = 0;
  else
    ConUpdateLinearCombination(linear_combination, constraints.back(), random_generator2.Uniform());

  /* Build the other linear constraints (if any) */
  for (j = 2, i = 1; j <= number_of_linear_constraints; ++j) {

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
        coco_error("CLinearConsBbobConstraints(): no inactive left i=%zu j=%zu nb_act=%zu nb_con=%zu", i, j,
                   number_of_active_constraints, number_of_linear_constraints);
      gradient = {};
      shift_factor = 0.01 + 2.0 * random_generator.Uniform();
    }

    constraints.push_back(CLinearSingleConstraint(
        function, dimension, instance,
        i <= number_of_active_constraints ? i : number_of_linear_constraints - inactive_constraints_left, factor1,
        gradient, shift_factor, feasible_direction));
    if (shift_factor == 0) {  /* active constraint */
      if (!gradient.empty()) { /* preserve the constraint based on function gradient */
        if (first_constraint != SIZE_MAX)
          coco_warning("CLinearConsBbobConstraints(): first_constraint already assigned, this is probably a bug");
        first_constraint = constraints.size() - 1;
      } else
        ConUpdateLinearCombination(linear_combination, constraints.back(), random_generator2.Uniform());
    }
  }

  /* Modify first constraint without changing the feasible solution,
   * thereby disguising the gradient of the function
   */
  if (number_of_active_constraints > 1 && disguise)
    constraints[first_constraint].DisguiseGradient(linear_combination, random_generator2.Uniform());

  return constraints;
}

} // namespace coco
