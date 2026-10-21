/**
 * @file suite_largescale.cpp
 * @brief Implementation of the bbob large-scale suite containing 24 functions in 6 large dimensions.
 */

#include "coco.h"

#include "suite_bbob_problems.cpp"

namespace coco {

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
  SuiteLargescale() : Suite("bbob-largescale", 24, {20, 40, 80, 160, 320, 640}, "instances: 1-15", true) {}

  std::string_view GetInstancesByYear(const int year) const override {
    if (year >= 2016 || year == 0) {
      return "1-15";
    } else if (year >= 2009) {
      return "1-15";
    } else {
      coco_error("SuiteLargescale::GetInstancesByYear(): year %d not defined for suite_largescale", year);
      return {};
    }
  }

  /** @brief Returns the problem with the given function, dimension and instance. */
  ProblemPtr GetProblem(const size_t function, const size_t dimension, const size_t instance) override {
    ProblemPtr problem;

    problem = GetLargescaleProblem(function, dimension, instance);

    return problem;
  }
};

} // namespace coco
