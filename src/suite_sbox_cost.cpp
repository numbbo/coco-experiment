/** @file suite_sbox_cost.cpp
 *
 * @brief Implementation of the bbob-boxed (aka SBOX-COST) suite, containing 
 * the 24 noiseless single-objective functions in 6 dimensions from the 
 * original BBOB suite.
 *
 * FIXME: Drop f9, f19, f20 and f24 which have optima that are always close to
 *   the origin?
 */

#include "coco.h"

#include "suite_bbob_problems.cpp"

#include "transform_vars_enforce_box.cpp"

namespace coco {

/**
 * @brief Sets the instances associated with years for the SBOX-COST suite.
 */
static std::string_view SuiteSboxCostGetInstancesByYear(const int year) {
  if (year == 2023) {
    return "1-5,101-110";
  } if (year == 0000) {
    return "1-15";
  } 
  else {
    coco_error(
        "suite_bbob_boxed_get_instances_by_year(): year %d not defined for "
        "suite_bbob_boxed",
        year);
    return {};
  }
}

/**
 * @brief Creates and returns a bbob-boxed problem without needing the actual
 * bbob-boxed suite.
 */
static ProblemPtr GetBbobBoxedProblem(const size_t function, const size_t dimension, const size_t instance) {
  ProblemPtr problem =
      MakeBbobProblem(function, dimension, instance, BbobVariant::BOXED, "bbob-boxed_f%03zu_i%02zu_d%02zu",
                      "bbob-boxed suite problem f%zu instance %zu in %zuD");
  return std::make_unique<TransformVarsEnforceBox>(std::move(problem));
}

/**
 * @brief The bbob-boxed suite.
 */
class SuiteSboxCost final : public Suite {
public:
  using Suite::Suite;

  std::string_view GetInstancesByYear(const int year) const override { return SuiteSboxCostGetInstancesByYear(year); }

  /** @brief Returns the problem with the given function, dimension and instance. */
  ProblemPtr GetProblem(const size_t function, const size_t dimension, const size_t instance) override {
    ProblemPtr problem;

    problem = GetBbobBoxedProblem(function, dimension, instance);

    return problem;
  }
};

/**
 * @brief Sets the dimensions and default instances for the bbob-boxed suite.
 */
static Suite *SuiteSboxCostInitialize() {

  Suite *suite;
  const size_t dimensions[] = {2, 3, 5, 10, 20, 40};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  suite = new SuiteSboxCost("bbob-boxed", 24, num_dimensions, dimensions, "year: 2023", 1);

  return suite;
}

} // namespace coco
