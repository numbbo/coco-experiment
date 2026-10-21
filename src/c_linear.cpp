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
static coco_problem_ptr_t c_guarantee_feasible_point(coco_problem_ptr_t problem, const double *feasible_point);

static coco_problem_ptr_t c_sum_variables_allocate(const size_t number_of_variables);

static coco_problem_ptr_t c_linear_transform(coco_problem_ptr_t inner_problem, const double *gradient,
                                             const double x_shift_factor);

static coco_problem_ptr_t
c_linear_single_cons_bbob_problem_allocate(const size_t function, const size_t dimension, const size_t instance,
                                           const size_t constraint_number, const double factor1,
                                           const char *problem_id_template, const char *problem_name_template,
                                           double *gradient, double x_shift_factor, const double *feasible_direction);

static coco_problem_ptr_t
c_linear_cons_bbob_problem_allocate(const size_t function, const size_t dimension, const size_t instance,
                                    const size_t number_of_active_constraints, const char *problem_id_template,
                                    const char *problem_name_template, const double *feasible_direction);

/**
 * @brief A linear constraint with all-ones gradient.
 */
struct c_sum_variables_t final : public coco_problem_s {
  explicit c_sum_variables_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 0, 1) {
    size_t i;
    set_id("%s_d%02lu", "linearconstraint", number_of_variables);
    for (i = 0; i < number_of_variables; ++i) {
      smallest_values_of_interest[i] = -5.0;
      largest_values_of_interest[i] = 5.0;
    }
  }

  /**
   * @brief Evaluates the linear constraint with all-ones gradient at
   *        the point 'x' and stores the result into 'y'.
   */
  void evaluate_constraint(const double *x, double *y, int update_counter) override {
    size_t i;

    assert(number_of_constraints == 1);
    y[0] = 0.0;
    for (i = 0; i < number_of_variables; ++i)
      /* prevent that the optimal solution is infeasible due to loss of precision when shifted */
      if (x[i] > 1e-11 || x[i] < 0)
        y[0] += x[i];

    (void)update_counter; /* To silence the compiler */
  }
};

/**
 * @brief A linear constraint whose gradient is given and whose feasible domain is extended by a shift factor.
 */
struct c_linear_t final : public coco_problem_transformed_t {
  std::vector<double> gradient;
  std::vector<double> transformed_x;
  const double x_shift_factor; /* shift solution by - factor * gradient */
  const double gradient_norm;

  c_linear_t(coco_problem_ptr_t inner, const double *gradient, const double x_shift_factor, const double gradient_norm)
      : coco_problem_transformed_t(std::move(inner), "gradient_linear_constraint"),
        gradient(gradient, gradient + inner_problem->number_of_variables),
        transformed_x(inner_problem->number_of_variables), x_shift_factor(x_shift_factor),
        gradient_norm(gradient_norm) {}

  /**
   * @brief Evaluates the linear constraint at the point 'x' and stores
   *        the result in 'y'.
   */
  void evaluate_constraint(const double *x, double *y, int update_counter) override {
    size_t i;
    double factor;

    assert(number_of_constraints == 1);
    for (i = 0, factor = x_shift_factor / gradient_norm; i < number_of_variables; ++i)
      transformed_x[i] = (gradient[i]) * (x[i] - factor * gradient[i]);

    inner_problem->evaluate_constraint(transformed_x.data(), y, update_counter);
  }
};

/**
 * @brief Returns the linear constraint behind the given problem pointer.
 */
static c_linear_t *c_linear_get(const coco_problem_t *problem) {
  assert(dynamic_cast<const c_linear_t *>(problem) != nullptr);
  return static_cast<c_linear_t *>(const_cast<coco_problem_t *>(problem));
}

/**
 * @brief Guarantees that "feasible_direction" is feasible w.r.t.
 *        the constraint in "problem" and records it as the
 *        initial feasible solution to this coco_problem.
 */
static coco_problem_ptr_t c_guarantee_feasible_point(coco_problem_ptr_t problem, const double *feasible_direction) {

  size_t i;
  c_linear_t *data;
  double constraint_value = 0.0;

  data = c_linear_get(problem.get());

  assert(problem->number_of_constraints == 1);

  /* Let p be the gradient of the constraint in "problem".
   * Check whether p' * (feasible_direction) <= 0.
   */
  problem->evaluate_constraint(feasible_direction, &constraint_value, 0);

  /* Flip the constraint in "problem" if feasible_direction
   * is not feasible w.r.t. the constraint in "problem".
   */
  if (constraint_value > 0)
    for (i = 0; i < problem->number_of_variables; ++i)
      data->gradient[i] *= -1.0;

  problem->initial_solution.assign(feasible_direction, feasible_direction + problem->number_of_variables);

  data = nullptr;
  return problem;
}

/**
 * @brief Allocates a linear constraint coco_problem_t with all-ones gradient.
 */
static coco_problem_ptr_t c_sum_variables_allocate(const size_t number_of_variables) {
  return std::make_unique<c_sum_variables_t>(number_of_variables);
}

/**
 * @brief Transforms a linear constraint with all-ones gradient
 *        into a linear constraint whose gradient is passed
 *        as argument and with the feasible domain extended by
 *        shift_factor.
 */
static coco_problem_ptr_t c_linear_transform(coco_problem_ptr_t inner_problem, const double *gradient,
                                             const double shift_factor) {

  double gradient_norm = coco_vector_norm(gradient, inner_problem->number_of_variables);

  if (gradient_norm <= 0)
    coco_error("c_linear_transform(): gradient norm %f<=0 zero", gradient_norm);

  return std::make_unique<c_linear_t>(std::move(inner_problem), gradient, shift_factor, gradient_norm);
}

/**
 * @brief Builds a coco_problem_t containing one single linear constraint.
 *
 * This function is called by c_linear_cons_bbob_problem_allocate(),
 * the central function that stacks all the constraints built by
 * c_linear_single_cons_bbob_problem_allocate() into one single
 * coco_problem_t object.
 */
static coco_problem_ptr_t c_linear_single_cons_bbob_problem_allocate(
    const size_t function, const size_t dimension, const size_t instance, const size_t constraint_number,
    const double factor1, const char *problem_id_template, const char *problem_name_template, double *gradient,
    const double x_shift_factor, const double *feasible_direction) {

  size_t i;

  coco_problem_ptr_t problem;
  long seed_cons_i;
  double factor2;

  problem = c_sum_variables_allocate(dimension);

  seed_cons_i = (long)(function + 10000 * instance + 50000 * constraint_number);
  coco_random_state_t random_generator((uint32_t)seed_cons_i);

  /* The constraints gradients are scaled with random numbers
   * 10**U[0,1] and 10**U_i[0,2], where U[a, b] is uniform in [a,b]
   * and only U_i is drawn for each constraint individually.
   * The random number 10**U[0,1] is given by the variable 'factor1'
   * while the random number 10**U_i[0,2] is calculated below and
   * stored as 'factor2'. (The exception is when the number of
   * constraints is n+1, in which case 'factor2' defines a random
   * number 10**U_i[0,1])
   */

  factor2 = pow(100.0, random_generator.uniform());

  /* Set the gradient of the linear constraint if it is given.
   * This should be the case of the construction of the first
   * linear constraint only.
   */
  if (gradient) {

    coco_vector_scale(gradient, dimension, factor1 * factor2, coco_vector_norm(gradient, dimension));
    problem = c_linear_transform(std::move(problem), gradient, x_shift_factor);

  } else { /* Randomly generate the gradient of the linear constraint */

    std::vector<double> gradient_linear_constraint(dimension);

    /* Generate a pseudorandom vector with distribution N_i(0, I)
     * and scale it with 'factor1' and 'factor2' (see comments above)
     */
    for (i = 0; i < dimension; ++i)
      gradient_linear_constraint[i] = factor1 * random_generator.normal() * factor2 / sqrt((double)dimension);

    problem = c_linear_transform(std::move(problem), gradient_linear_constraint.data(), x_shift_factor);
  }

  /* Guarantee that the vector feasible_point is feasible w.r.t. to
   * this constraint and set it as the initial solution.
   * The initial solution will be copied later to the constrained function
   * coco_problem_t object once the objective function and the constraint(s)
   * are stacked together in coco_problem_stacked_allocate().
   */
  if (feasible_direction)
    problem = c_guarantee_feasible_point(std::move(problem), feasible_direction);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("linear");
  return problem;
}

/**
 * @brief helper function to successively compute a linear combination of
 *        constraint gradients, add current gradient if weight is != 0.
 *
 */
static void con_update_linear_combination(double *linear_combination, const coco_problem_t *problem, double weight) {
  size_t i;
  const c_linear_t *data;

  data = c_linear_get(problem);
  if (data->gradient.empty()) {
    if (weight != 0) {
      coco_error("con_update_linear_combination(): gradient of constraint was zero");
    } else {
      coco_warning("con_update_linear_combination(): gradient of constraint was zero");
    }
  } else {
    if (data->x_shift_factor != 0)
      coco_warning("Inactive constraint passed to update_linear_combination, x_shift_factor=%f", data->x_shift_factor);
    if (weight == 0)
      return; /* nothing to add */
    if (weight < 0)
      coco_warning("con_update_linear_combination: weight=%f < 0, should be > 0", weight);
    for (i = 0; i < problem->number_of_variables; ++i)
      linear_combination[i] += weight * data->gradient[i];
  }
}

/**
 * @brief Builds a coco_problem_t containing all the linear constraints
 *        by stacking them all.
 *
 * The constraints' gradients are randomly generated with distribution
 * 10**U[0,1] * N_i(0, I) * 10**U_i[0,2], where U[a, b] is uniform
 * in [a,b] and only U_i is drawn for each constraint individually.
 * The exception is the first constraint, whose gradient is given by
 * 10**U[0,1] * (-feasible_direction) * 10**U_i[0,2].
 *
 * Each constraint is built by calling the function
 * c_linear_single_cons_bbob_problem_allocate(), which returns a
 * coco_problem_t object that defines the constraint. The resulting
 * coco_problem_t objects are then stacked together into one single
 * coco_problem_t object that is returned by the function.
 */
static coco_problem_ptr_t
c_linear_cons_bbob_problem_allocate(const size_t function, const size_t dimension, const size_t instance,
                                    const size_t number_of_active_constraints, const char *problem_id_template,
                                    const char *problem_name_template, const double *feasible_direction) {

  const double global_scaling_factor = 100.;
  size_t i, j;

  coco_problem_ptr_t problem_c;
  coco_problem_ptr_t problem_c2;
  double *gradient;
  c_linear_t *first_constraint_data = nullptr;
  double shift_factor;
  long seed_cons;
  double fac, norm, factor1;
  size_t number_of_linear_constraints;
  size_t inactive_constraints_left =
      number_of_active_constraints / 2; /* TODO: decide whether this should go into the interface */
  int disguise_gradient = 1;

  if (std::string_view(problem_id_template).substr(0, 28) == "bbob-constrained-active-only")
    inactive_constraints_left = 0;

  if (std::string_view(problem_id_template).substr(0, 28) == "bbob-constrained-no-disguise")
    disguise_gradient = 0;

  number_of_linear_constraints = number_of_active_constraints + inactive_constraints_left;

  std::vector<double> gradient_c1(dimension);
  std::vector<double> linear_combination(dimension, 0.0);

  for (i = 0; i < dimension; ++i)
    gradient_c1[i] = -feasible_direction[i];

  /* Build a coco_problem_t object for each constraint.
   * The constraints' gradients are generated randomly with
   * distribution 10**U[0,1] * N_i(0, I/n) * 10**U_i[0,2]
   * where U[a, b] is uniform in [a,b] and only U_i is drawn
   * for each constraint individually.
   */

  /* Calculate the first random factor 10**U[0,1]. */
  seed_cons = (long)(function + 10000 * instance);
  coco_random_state_t random_generator((uint32_t)seed_cons);
  coco_random_state_t random_generator2((uint32_t)seed_cons);
  factor1 = global_scaling_factor * pow(10.0, random_generator.uniform());

  /* Build the first linear constraint using 'gradient_c1' to build
   * its gradient.
   */
  /* set gradient depending on instance number */
  gradient = instance % number_of_active_constraints ? nullptr : gradient_c1.data();
  problem_c = c_linear_single_cons_bbob_problem_allocate(function, dimension, instance, 1, factor1, problem_id_template,
                                                         problem_name_template, gradient, 0.0, feasible_direction);
  if (gradient != nullptr) /* preserve the constraint based on function gradient */
    first_constraint_data = c_linear_get(problem_c.get());
  else
    con_update_linear_combination(linear_combination.data(), problem_c.get(), random_generator2.uniform());

  /* Instantiate the other linear constraints (if any) and stack them
   * all into problem_c
   */
  for (j = 2, i = 1; j <= number_of_linear_constraints; ++j) {

    /* Instantiate a new problem containing one linear constraint only */
    /* set gradient depending on instance number */

    if (i < number_of_active_constraints &&
        random_generator.uniform() + 1e-23 >
            (double)inactive_constraints_left /
                (double)number_of_linear_constraints) { /* create an active constraint */
      gradient = (++i - 1 + instance) % number_of_active_constraints ? nullptr : gradient_c1.data();
      shift_factor = 0.0;
    } else { /* create an inactive (shifted) constraint */
      if (inactive_constraints_left-- <= 0)
        coco_error("c_linear_cons_bbob_problem_allocate(): no inactive left i=%lu j=%lu nb_act=%lu nb_con=%lu",
                   (unsigned long)i, (unsigned long)j, (unsigned long)number_of_active_constraints,
                   (unsigned long)number_of_linear_constraints);
      gradient = nullptr;
      shift_factor = 0.01 + 2.0 * random_generator.uniform();
    }

    problem_c2 = c_linear_single_cons_bbob_problem_allocate(
        function, dimension, instance,
        i <= number_of_active_constraints ? i : number_of_linear_constraints - inactive_constraints_left, factor1,
        problem_id_template, problem_name_template, gradient, shift_factor, feasible_direction);
    if (shift_factor == 0) {  /* active constraint */
      if (gradient != nullptr) { /* preserve the constraint based on function gradient */
        if (first_constraint_data)
          coco_warning(
              "c_linear_cons_bbob_problem_allocate(): first_constraint_data already assigned, this is probably a bug");
        first_constraint_data = c_linear_get(problem_c2.get());
      } else
        con_update_linear_combination(linear_combination.data(), problem_c2.get(), random_generator2.uniform());
    }

    const std::string problem_type_c2 = problem_c2->problem_type;
    problem_c = coco_problem_stacked_allocate(std::move(problem_c), std::move(problem_c2));

    /* Use the standard stacked problem_id as problem_name and
     * construct a new suite-specific problem_id
     */
    problem_c->set_name("%s", problem_c->problem_id.c_str());
    problem_c->set_id("bbob-constrained_f%02lu_i%02lu_d%02lu", (unsigned long)function, (unsigned long)instance,
                      (unsigned long)dimension);

    /* Construct problem type */
    problem_c->set_type("%s_%s", problem_type_c2.c_str(), problem_type_c2.c_str());
  }

  /* Modify first constraint without changing the feasible solution,
   * thereby disguising the gradient of the function
   */
  if (number_of_active_constraints > 1 && disguise_gradient) {
    norm = first_constraint_data->gradient_norm; /* for reading convenience only */
    gradient = first_constraint_data->gradient.data(); /* ditto */
    fac = coco_vector_scalar_product(linear_combination.data(), gradient, dimension);
    if (fac < 0)
      coco_error("scalar product between first constraint and linear combination = %f < 0 should be > 0", fac);
    fac = fac <= 0 ? (double)dimension : coco_double_min((double)dimension, norm * norm / fac);
    fac *= (1 + random_generator2.uniform()) / (2 + 1. / (double)dimension); /* bound to <1 and >1/3 */
    for (i = 0; i < dimension; ++i)
      gradient[i] -= fac * linear_combination[i];
    fac = coco_vector_norm(gradient, dimension); /* new gradient norm */
    if (fac == 0)                                /* make sure that gradient did not become zero */
      coco_error("new vector norm = %f == 0, should be > 0", fac);
    coco_vector_scale(gradient, dimension, norm, fac); /* restore original norm */
  }


  return problem_c;
}
