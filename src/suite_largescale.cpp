/**
 * @file suite_largescale.cpp
 * @brief Implementation of the bbob large-scale suite containing 24 functions in 6 large dimensions.
 */

#include "coco.h"

#include "suite_bbob_problems.cpp"

namespace coco {

/**
 * @brief Sets the instances associated with years for the bbob large-scale suite.
 */
static std::string_view SuiteLargescaleGetInstancesByYear(const int year) {
  if (year >= 2016 || year == 0) {
    return "1-15";
  } else if (year >= 2009) {
    return "1-15";
  } else {
    coco_error("SuiteLargescaleGetInstancesByYear(): year %d not defined for suite_largescale", year);
    return {};
  }
}

/**
 * @brief Creates and returns a large-scale problem without needing the actual large-scale suite.
 */
static ProblemPtr GetLargescaleProblem(const size_t function, const size_t dimension, const size_t instance) {
  return MakeBbobProblem(function, dimension, instance, BbobVariant::LARGESCALE, "bbob_f%03zu_i%02zu_d%04zu",
                         "BBOB large-scale suite problem f%zu instance %zu in %zuD");
}

/**
 * @brief The bbob-largescale suite.
 */
class SuiteLargescale final : public Suite {
public:
  using Suite::Suite;

  std::string_view GetInstancesByYear(const int year) const override { return SuiteLargescaleGetInstancesByYear(year); }

  /** @brief Returns the problem with the given function, dimension and instance. */
  ProblemPtr GetProblem(const size_t function, const size_t dimension, const size_t instance) override {
    ProblemPtr problem;

    problem = GetLargescaleProblem(function, dimension, instance);

    return problem;
  }
};

/**
 * @brief Sets the dimensions and default instances for the bbob large-scale suite.
 */
static Suite *SuiteLargescaleInitialize() {

  Suite *suite;
  const size_t dimensions[] = {20, 40, 80, 160, 320, 640};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  suite = new SuiteLargescale("bbob-largescale", 24, num_dimensions, dimensions, "instances: 1-15", 1);
  return suite;
}

} // namespace coco
