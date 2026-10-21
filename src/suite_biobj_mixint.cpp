/**
 * @file suite_biobj_mixint.cpp
 * @brief Implementation of a bi-objective mixed-integer suite. The functions are the same as those
 * in the bbob-biobj-ext suite with 92 functions, but the large-scale implementations of the
 * functions are used instead of the original ones for dimensions over 40.
 */

#include "coco.h"
#include "mo_utilities.cpp"
#include "suite_biobj_utilities.cpp"
#include "suite_largescale.cpp"
#include "transform_vars_discretize.cpp"
#include "transform_obj_scale.cpp"
#include "suite_bbob_mixint.cpp"

/**
 * @brief Sets the instances associated with years for the bbob-biobj-mixint suites.
 *
 * @note The instances of the bi-objective suites generally do not changes with years.
 */
static std::string_view suite_biobj_mixint_get_instances_by_year(const int year) {

  (void)year; /* To get rid of compiler warnings */
  return "1-15";
}

/**
 * @brief Creates and returns a mixed-integer bi-objective bbob problem without needing the actual
 * bbob-mixint suite.
 *
 * The problem is constructed by first finding the underlying single-objective continuous problems,
 * then discretizing the problems, then scaling them to adjust their difficulty and finally stacking
 * them to get a bi-objective mixed-integer problem.
 *
 * @param function Function
 * @param dimension Dimension
 * @param instance Instance
 * @param coco_get_problem_function The function that is used to access the single-objective problem.
 * @param new_instances The new instances computed so far (new ones are added).
 * @param dimensions An array of dimensions to take into account when creating new instances.
 * @param num_dimensions The number of dimensions to take into account when creating new instances.
 * @return The problem that corresponds to the given parameters.
 */
static coco_problem_ptr_t coco_get_biobj_mixint_problem(const size_t function, const size_t dimension,
                                                        const size_t instance,
                                                        const coco_get_problem_function_t coco_get_problem_function,
                                                        suite_biobj_new_instances_t &new_instances,
                                                        const size_t *dimensions, const size_t num_dimensions) {

  coco_problem_ptr_t problem_cont;

  coco_problem_ptr_t problem;
  coco_problem_ptr_t problem1;
  coco_problem_ptr_t problem2;
  coco_problem_ptr_t problem1_mixint;
  coco_problem_ptr_t problem2_mixint;

  std::vector<double> smallest_values_of_interest(dimension);
  std::vector<double> largest_values_of_interest(dimension);

  size_t j;
  size_t num_integer = dimension;
  /* The cardinality of variables (0 = continuous variables should always come last) */
  const size_t variable_cardinality[] = {2, 4, 8, 16, 0};

  if (dimension % 5 != 0)
    coco_error("coco_get_biobj_mixint_problem(): dimension %lu not supported for suite_bbob_mixint", dimension);

  /* First, find the underlying single-objective continuous problems */
  problem_cont = coco_get_biobj_problem(function, dimension, instance, coco_get_problem_function, new_instances,
                                        dimensions, num_dimensions);
  assert(problem_cont != nullptr);
  /* Take the two single-objective problems out of the stacked problem and free only the latter */
  problem1 = std::move(coco_problem_as_stacked(problem_cont.get())->problem1);
  problem2 = std::move(coco_problem_as_stacked(problem_cont.get())->problem2);
  assert(problem1);
  assert(problem2);
  problem_cont.reset();

  /* Set the ROI of the outer problem according to the given cardinality of variables and the ROI of the
   * inner problems to [-4, 4] for variables that will be discretized */
  for (size_t i = 0; i < dimension; i++) {
    j = i / (dimension / 5);
    if (variable_cardinality[j] == 0) {
      /* Continuous variables */
      /* Outer problem */
      smallest_values_of_interest[i] = -100;
      largest_values_of_interest[i] = 100;
      if (num_integer == dimension)
        num_integer = i;
    } else {
      /* Outer problem */
      smallest_values_of_interest[i] = 0;
      largest_values_of_interest[i] = (double)variable_cardinality[j] - 1;
      /* Inner problems */
      problem1->smallest_values_of_interest[i] = -4;
      problem1->largest_values_of_interest[i] = 4;
      problem2->smallest_values_of_interest[i] = -4;
      problem2->largest_values_of_interest[i] = 4;
    }
  }

  const size_t function1 = coco_problem_get_suite_dep_function(problem1.get());
  const size_t function2 = coco_problem_get_suite_dep_function(problem2.get());
  const std::string problem_type1 = problem1->problem_type;
  const std::string problem_type2 = problem2->problem_type;

  /* Second, discretize the single-objective problems */
  problem1_mixint = transform_vars_discretize(std::move(problem1), smallest_values_of_interest.data(),
                                              largest_values_of_interest.data(), num_integer);
  problem2_mixint = transform_vars_discretize(std::move(problem2), smallest_values_of_interest.data(),
                                              largest_values_of_interest.data(), num_integer);

  /* Third, scale the objective values */
  problem1_mixint = transform_obj_scale(std::move(problem1_mixint), suite_bbob_mixint_scaling_factors[function1 - 1]);
  problem2_mixint = transform_obj_scale(std::move(problem2_mixint), suite_bbob_mixint_scaling_factors[function2 - 1]);

  /* Fourth, combine the problems in a bi-objective mixed-integer problem */
  problem = coco_problem_stacked_allocate(std::move(problem1_mixint), std::move(problem2_mixint),
                                          smallest_values_of_interest.data(), largest_values_of_interest.data());

  /* Use the standard stacked problem_id as problem_name and construct a new problem_id */
  problem->set_name("%s", problem->problem_id.c_str());
  problem->set_id("bbob-biobj-mixint_f%03lu_i%02lu_d%03lu", (unsigned long)function, (unsigned long)instance,
                  (unsigned long)dimension);

  /* Construct problem type */
  problem->set_type("%s_%s", problem_type1.c_str(), problem_type2.c_str());

  return problem;
}
/**
 * @brief The bbob-biobj-mixint suite.
 */
struct suite_biobj_mixint_t final : public coco_suite_s {
  using coco_suite_s::coco_suite_s;

  std::string_view get_instances_by_year(const int year) const override {
    return suite_biobj_mixint_get_instances_by_year(year);
  }

  /** @brief Returns the problem with the given function, dimension and instance. */
  coco_problem_ptr_t get_problem(const size_t function, const size_t dimension, const size_t instance) override {
    coco_problem_ptr_t problem;
    const size_t dim_large_scale = 50; /* Switch to large-scale functions for dimensions over 50 */

    if (dimension < dim_large_scale)
      problem = coco_get_biobj_mixint_problem(function, dimension, instance, coco_make_bbob_problem, new_instances,
                                              dimensions.data(), dimensions.size());
    else
      problem = coco_get_biobj_mixint_problem(function, dimension, instance, coco_get_largescale_problem, new_instances,
                                              dimensions.data(), dimensions.size());

    return problem;
  }

  suite_biobj_new_instances_t new_instances; /**< @brief The new instances computed so far. */
};

/**
 * @brief Sets the dimensions and default instances for the bbob-biobj-mixint suite.
 */
static coco_suite_t *suite_biobj_mixint_initialize(void) {

  coco_suite_t *suite;
  const size_t dimensions[] = {5, 10, 20, 40, 80, 160};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  suite = new suite_biobj_mixint_t("bbob-biobj-mixint", 92, num_dimensions, dimensions, "instances: 1-15", 0);

  return suite;
}
