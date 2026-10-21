/**
 * @file  suite_cons_bbob.cpp
 * @brief Implementation of the constrained bbob suite containing
 *        54 constrained problems in 6 dimensions. See comments in
 *        "suite_cons_bbob_problems.cpp" for more details.
 */

#include "coco.h"
#include "suite_cons_bbob_problems.cpp"

namespace coco {

/**
 * @brief Sets the instances associated with years for the constrained
 *        bbob suite.
 */
static std::string_view SuiteConsBbobGetInstancesByYear(const int year) {

  if ((year >= 2022) || (year == 0)) {
    return "1-15";
  } else if ((year >= 2009)) {
    return "1-15";
  } else {
    coco_error("SuiteConsBbobGetInstancesByYear(): year %d not defined for suite_cons_bbob", year);
    return {};
  }
}

/**
 * @brief The bbob-constrained suites.
 */
class SuiteConsBbob final : public Suite {
public:
  using Suite::Suite;

  std::string_view GetInstancesByYear(const int year) const override { return SuiteConsBbobGetInstancesByYear(year); }

  /** @brief Returns the problem with the given function, dimension and instance. */
  ProblemPtr GetProblem(const size_t function, const size_t dimension, const size_t instance) override {
    ProblemPtr problem;

    problem = GetConsBbobProblem(name(), function, dimension, instance);

    /* Construct a suite-specific problem id and name */
    problem->SetId("%s_f%03zu_i%02zu_d%02zu", name().c_str(), function, instance, dimension);
    problem->SetName("%s suite problem f%zu instance %zu in %zuD", name().c_str(), function, instance, dimension);

    return problem;
  }
};

/**
 * @brief Sets the dimensions and default instances for the bbob suite.
 */
static Suite *SuiteConsBbobInitialize(std::string_view suite_name) {

  Suite *suite;
  const size_t dimensions[] = {2, 3, 5, 10, 20, 40};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  /* IMPORTANT: Make sure to change the default instance for every new workshop! */
  suite = new SuiteConsBbob(suite_name, 54, num_dimensions, dimensions, "year: 2022", 1);

  return suite;
}

} // namespace coco
