/**
 * @file suite_biobj.cpp
 * @brief Implementation of two bi-objective suites created by combining two single-objective problems
 * from the bbob suite:
 * - bbob-biobj contains 55 functions and 6 dimensions
 * - bbob-biobj-ext contains 55 + 37 functions and 6 dimensions
 *
 * The 55 functions of the bbob-biobj suite are created by combining any two single-objective bbob functions
 * i,j (where i<j) from a subset of 10 functions.
 *
 * The first 55 functions of the bbob-biobj-ext suite are the same as in the original bbob-biobj test suite
 * to which 37 functions are added. Those additional functions are constructed by combining all not yet
 * contained in-group combinations (i,j) of single-objective bbob functions i and j such that i<j (i.e. in
 * particular not all combinations (i,i) are included in this bbob-biobj-ext suite), with the exception of
 * the Weierstrass function (f16) for which the optimum is not unique and thus a nadir point is difficult
 * to compute, see http://numbbo.github.io/coco-doc/bbob-biobj/functions/ for details.
 *
 * @note See file suite_biobj_utilities.cpp for the implementation of the bi-objective problems and the handling
 * of new instances.
 */

#include "coco.h"
#include "mo_utilities.cpp"
#include "suite_biobj_utilities.cpp"
#include "suite_bbob.cpp"

/**
 * @brief Sets the instances associated with years for the bbob-biobj suites.
 *
 * @note The instances of the bi-objective suites generally do not changes with years.
 */
static std::string_view suite_biobj_get_instances_by_year(const int year) {

  if ((year == 2016) || (year == 0000)) { /* test case */
    return "1-10";
  } else
    return "1-15";
}

/**
 * @brief The bbob-biobj and bbob-biobj-ext suites.
 */
struct suite_biobj_t final : public coco_suite_s {
  using coco_suite_s::coco_suite_s;

  std::string_view get_instances_by_year(const int year) const override {
    return suite_biobj_get_instances_by_year(year);
  }

  /** @brief Returns the problem with the given function, dimension and instance. */
  coco_problem_ptr_t get_problem(const size_t function, const size_t dimension, const size_t instance) override {
    coco_problem_ptr_t problem;

    problem = coco_get_biobj_problem(function, dimension, instance, coco_make_bbob_problem, new_instances,
                                     dimensions.data(), dimensions.size());

    problem->is_opt_known = 1;

    return problem;
  }

  suite_biobj_new_instances_t new_instances; /**< @brief The new instances computed so far. */
};

/**
 * @brief Sets the dimensions and default instances for the bbob-biobj suites.
 */
static coco_suite_t *suite_biobj_initialize(std::string_view suite_name) {

  coco_suite_t *suite;
  const size_t dimensions[] = {2, 3, 5, 10, 20, 40};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  if (suite_name == "bbob-biobj") {
    suite = new suite_biobj_t("bbob-biobj", 55, num_dimensions, dimensions, "instances: 1-15", 1);
  } else if (suite_name == "bbob-biobj-ext") {
    suite = new suite_biobj_t("bbob-biobj-ext", 55 + 37, num_dimensions, dimensions, "instances: 1-15", 1);
  } else {
    coco_error("suite_biobj_initialize(): unknown problem suite");
    return nullptr;
  }

  return suite;
}
