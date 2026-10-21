/**
 * @file suite_largescale.cpp
 * @brief Implementation of the bbob large-scale suite containing 24 functions in 6 large dimensions.
 */

#include "coco.h"

#include "suite_bbob_problems.cpp"

/**
 * @brief Sets the instances associated with years for the bbob large-scale suite.
 */
static std::string_view suite_largescale_get_instances_by_year(const int year) {
  if (year >= 2016 || year == 0) {
    return "1-15";
  } else if (year >= 2009) {
    return "1-15";
  } else {
    coco_error("suite_largescale_get_instances_by_year(): year %d not defined for suite_largescale", year);
    return {};
  }
}

/**
 * @brief Creates and returns a large-scale problem without needing the actual large-scale suite.
 */
static coco_problem_ptr_t coco_get_largescale_problem(const size_t function, const size_t dimension,
                                                      const size_t instance) {
  return bbob_make_problem(function, dimension, instance, coco_bbob_variant_t::largescale, "bbob_f%03zu_i%02zu_d%04zu",
                           "BBOB large-scale suite problem f%zu instance %zu in %zuD");
}

/**
 * @brief The bbob-largescale suite.
 */
class suite_largescale_t final : public coco_suite_s {
public:
  using coco_suite_s::coco_suite_s;

  std::string_view get_instances_by_year(const int year) const override {
    return suite_largescale_get_instances_by_year(year);
  }

  /** @brief Returns the problem with the given function, dimension and instance. */
  coco_problem_ptr_t get_problem(const size_t function, const size_t dimension, const size_t instance) override {
    coco_problem_ptr_t problem;

    problem = coco_get_largescale_problem(function, dimension, instance);

    return problem;
  }
};

/**
 * @brief Sets the dimensions and default instances for the bbob large-scale suite.
 */
static coco_suite_t *suite_largescale_initialize() {

  coco_suite_t *suite;
  const size_t dimensions[] = {20, 40, 80, 160, 320, 640};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  suite = new suite_largescale_t("bbob-largescale", 24, num_dimensions, dimensions, "instances: 1-15", 1);
  return suite;
}
