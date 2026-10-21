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
 * @brief The bbob-constrained suites.
 */
class SuiteConsBbob final : public Suite {
public:
  /** @brief Creates the constrained suite with the given name (bbob-constrained or one of its variants). */
  /* IMPORTANT: Make sure to change the default instance for every new workshop! */
  explicit SuiteConsBbob(std::string_view suite_name)
      : Suite(suite_name, 54, {2, 3, 5, 10, 20, 40}, "year: 2022", true) {}

  std::string_view GetInstancesByYear(const int year) const override {
    if ((year >= 2022) || (year == 0)) {
      return "1-15";
    } else if ((year >= 2009)) {
      return "1-15";
    } else {
      coco_error("SuiteConsBbob::GetInstancesByYear(): year %d not defined for suite_cons_bbob", year);
      return {};
    }
  }

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

} // namespace coco
