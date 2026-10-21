/**
 * @file suite_bbob.cpp
 * @brief Implementation of the bbob suite containing 24 noiseless single-objective functions in 6
 * dimensions.
 */

#include "coco.h"

#include "suite_bbob_problems.cpp"

namespace coco {

/**
 * @brief Sets the instances associated with years for the bbob suite.
 */
static std::string_view SuiteBbobGetInstancesByYear(const int year) {

  if (year >= 2023) {
    return "1-5,101-110";
  } else if (year >= 2021) {
    return "1-5,91-100";
  } else if (year >= 2018) {
    return "1-5,71-80";
  } else if (year == 2017) {
    return "1-5,61-70";
  } else if ((year == 2016) || (year == 0000)) { /* test case */
    return "1-5,51-60";
  } else if (year == 2015) {
    return "1-5,41-50";
  } else if (year >= 2013) {
    return "1-5,31-40";
  } else if (year == 2012) {
    return "1-5,21-30";
  } else if (year >= 2010) {
    return "1-15";
  } else if (year == 2009) {
    return "1-5,1-5,1-5";
  } else {
    coco_error("SuiteBbobGetInstancesByYear(): year %d not defined for suite_bbob", year);
    return {};
  }
}

/**
 * @brief Creates and returns a BBOB problem without needing the actual bbob suite.
 *
 * Useful for other suites as well (see for example suite_biobj.cpp).
 */
static ProblemPtr GetBbobProblem(const size_t function, const size_t dimension, const size_t instance) {
  return MakeBbobProblem(function, dimension, instance, BbobVariant::BBOB, "bbob_f%03zu_i%02zu_d%02zu",
                         "BBOB suite problem f%zu instance %zu in %zuD");
}

/**
 * @brief The bbob suite.
 */
class SuiteBbob final : public Suite {
public:
  using Suite::Suite;

  std::string_view GetInstancesByYear(const int year) const override { return SuiteBbobGetInstancesByYear(year); }

  /** @brief Returns the problem with the given function, dimension and instance. */
  ProblemPtr GetProblem(const size_t function, const size_t dimension, const size_t instance) override {
    ProblemPtr problem;

    problem = GetBbobProblem(function, dimension, instance);

    return problem;
  }
};

/**
 * @brief Sets the dimensions and default instances for the bbob suite.
 */
static Suite *SuiteBbobInitialize() {

  Suite *suite;
  const size_t dimensions[] = {2, 3, 5, 10, 20, 40};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  /* IMPORTANT: Make sure to change the default instance for every new workshop! */
  suite = new SuiteBbob("bbob", 24, num_dimensions, dimensions, "year: 2018", 1);

  return suite;
}

} // namespace coco

/* C interface (declared in coco.h) */

/**
 * @brief Creates and returns a BBOB problem, which the caller must free with coco_problem_free().
 *
 * Used by the language bindings.
 */
coco_problem_t *coco_get_bbob_problem(const size_t function, const size_t dimension, const size_t instance) {
  return coco::GetBbobProblem(function, dimension, instance).release();
}
