/**
 * @file suite_toy.cpp
 * @brief Implementation of a toy suite containing 6 noiseless "basic" single-objective functions in 5
 * dimensions.
 */

#include "coco_problem.cpp"
#include "f_bueche_rastrigin.cpp"
#include "f_ellipsoid.cpp"
#include "f_linear_slope.cpp"
#include "f_rastrigin.cpp"
#include "f_rosenbrock.cpp"
#include "f_sphere.cpp"

/**
 * @brief Returns the problem from the toy suite that corresponds to the given parameters.
 *
 * @param suite The COCO suite.
 * @param function_idx Index of the function (starting from 0).
 * @param dimension_idx Index of the dimension (starting from 0).
 * @param instance_idx Index of the instance (starting from 0).
 * @return The problem that corresponds to the given parameters.
 */
static coco_problem_ptr_t suite_toy_get_problem(coco_suite_t *suite, const size_t function_idx,
                                                const size_t dimension_idx, const size_t instance_idx) {

  coco_problem_ptr_t problem;

  const size_t function = suite->functions[function_idx];
  const size_t dimension = suite->dimensions[dimension_idx];
  const size_t instance = suite->instances[instance_idx];

  if (function == 1) {
    problem = f_sphere_allocate(dimension);
  } else if (function == 2) {
    problem = f_ellipsoid_allocate(dimension, 1.0e6);
  } else if (function == 3) {
    problem = f_rastrigin_allocate(dimension);
  } else if (function == 4) {
    problem = f_bueche_rastrigin_allocate(dimension);
  } else if (function == 5) {
    double xopt[40] = {5.0};
    problem = f_linear_slope_allocate(dimension, xopt);
  } else if (function == 6) {
    problem = f_rosenbrock_allocate(dimension);
  } else {
    coco_error("suite_toy_get_problem(): function %lu does not exist in this suite", (unsigned long)function);
    return nullptr; /* Never reached */
  }

  problem->suite_dep_function = function;
  problem->suite_dep_instance = instance;
  problem->suite_dep_index = coco_suite_encode_problem_index(suite, function_idx, dimension_idx, instance_idx);

  return problem;
}

/**
 * @brief The toy suite.
 */
struct suite_toy_t final : public coco_suite_s {
  using coco_suite_s::coco_suite_s;

  coco_problem_ptr_t get_problem(const size_t function_idx, const size_t dimension_idx,
                                 const size_t instance_idx) override {
    return suite_toy_get_problem(this, function_idx, dimension_idx, instance_idx);
  }
};

/**
 * @brief Sets the dimensions and default instances for the toy suite.
 */
static coco_suite_t *suite_toy_initialize(void) {

  coco_suite_t *suite;
  const size_t dimensions[] = {2, 3, 5, 10, 20};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  suite = new suite_toy_t("toy", 6, num_dimensions, dimensions, "instances: 1", 1);

  return suite;
}
