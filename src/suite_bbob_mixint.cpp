/**
 * @file suite_bbob_mixint.cpp
 * @brief Implementation of a suite with mixed-integer bbob problems. The functions are the same
 * as those in the bbob suite with 24 functions, but the large-scale implementations of the
 * functions are used instead of the original ones for dimensions over 40. Additionally, the
 * functions are scaled as given by suite_bbob_mixint_scaling_factors.
 */

#include "coco.h"
#include "suite_bbob.cpp"
#include "suite_largescale.cpp"
#include "transform_vars_discretize.cpp"
#include "transform_obj_scale.cpp"

/**
 * @brief Factors used to scale separate functions in order to achieve similar difficulty.
 */
static double suite_bbob_mixint_scaling_factors[] = {
    1,    1e-3, 1e-1, 1e-1, 1,    /* f1 to f5 */
    1e-2, 1,    1e-2, 1e-2, 1e-3, /* f6 to f10 */
    1e-2, 1e-4, 1e-1, 1,    1e-1, /* f11 to f15 */
    1,    10,   1,    10,   1e-1, /* f16 to f20 */
    1,    1,    10,   1e-1        /* f21 to f24 */
};

/**
 * @brief Sets the instances associated with years for the bbob-mixint suites.
 */
static std::string_view suite_bbob_mixint_get_instances_by_year(const int) { return "1-15"; }

/**
 * @brief Creates and returns a mixed-integer bbob problem without needing the actual bbob-mixint
 * suite.
 *
 * @param function Function
 * @param dimension Dimension
 * @param instance Instance
 * @param coco_get_problem_function The function that is used to access the continuous problem.
 * @return The problem that corresponds to the given parameters.
 */
static coco_problem_ptr_t coco_get_bbob_mixint_problem(const size_t function, const size_t dimension,
                                                       const size_t instance,
                                                       const coco_get_problem_function_t coco_get_problem_function) {
  coco_problem_ptr_t problem;

  /* The cardinality of variables (0 = continuous variables should always come last) */
  const size_t variable_cardinality[] = {2, 4, 8, 16, 0};

  std::vector<double> smallest_values_of_interest(dimension);
  std::vector<double> largest_values_of_interest(dimension);
  std::string inner_problem_id;

  size_t j;
  size_t cardinality = 0;
  size_t num_integer = dimension;
  if (dimension % 5 != 0)
    coco_error("coco_get_bbob_mixint_problem(): dimension %zu not supported for suite_bbob_mixint", dimension);

  problem = coco_get_problem_function(function, dimension, instance);
  assert(problem != nullptr);

  /* Set the ROI of the outer problem according to the given cardinality of variables and the ROI of the
   * inner problem to [-4, 4] for variables that will be discretized */
  for (size_t i = 0; i < dimension; i++) {
    j = i / (dimension / 5);
    cardinality = variable_cardinality[j];
    if (cardinality == 0) {
      /* Continuous variables */
      /* Outer problem */
      smallest_values_of_interest[i] = -5;
      largest_values_of_interest[i] = 5;
      if (num_integer == dimension)
        num_integer = i;
    } else {
      /* Outer problem */
      smallest_values_of_interest[i] = 0;
      largest_values_of_interest[i] = static_cast<double>(cardinality) - 1;
      /* Inner problem */
      problem->smallest_values_of_interest[i] = -4;
      problem->largest_values_of_interest[i] = 4;
    }
  }

  inner_problem_id = problem->problem_id;

  problem = transform_vars_discretize(std::move(problem), smallest_values_of_interest.data(),
                                      largest_values_of_interest.data(), num_integer);

  problem = transform_obj_scale(std::move(problem), suite_bbob_mixint_scaling_factors[function - 1]);

  problem->set_id("bbob-mixint_f%03zu_i%02zu_d%02zu", function, instance, dimension);
  problem->set_name("mixint(%s)", inner_problem_id.c_str());

  return problem;
}

/**
 * @brief The bbob-mixint suite.
 */
struct suite_bbob_mixint_t final : public coco_suite_s {
  using coco_suite_s::coco_suite_s;

  std::string_view get_instances_by_year(const int year) const override {
    return suite_bbob_mixint_get_instances_by_year(year);
  }

  /** @brief Returns the problem with the given function, dimension and instance. */
  coco_problem_ptr_t get_problem(const size_t function, const size_t dimension, const size_t instance) override {
    coco_problem_ptr_t problem;
    const size_t dim_large_scale = 50; /* Switch to large-scale functions for dimensions over 50 */

    if (dimension < dim_large_scale)
      problem = coco_get_bbob_mixint_problem(function, dimension, instance, coco_make_bbob_problem);
    else
      problem = coco_get_bbob_mixint_problem(function, dimension, instance, coco_get_largescale_problem);

    return problem;
  }
};

/**
 * @brief Sets the dimensions and default instances for the bbob-mixint suite.
 */
static coco_suite_t *suite_bbob_mixint_initialize(std::string_view suite_name) {

  coco_suite_t *suite;
  const size_t dimensions[] = {5, 10, 20, 40, 80, 160};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);
  suite = new suite_bbob_mixint_t(suite_name, 24, num_dimensions, dimensions, "instances: 1-15", 1);

  return suite;
}
