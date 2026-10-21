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

/**
 * @brief Sets the instances associated with years for the SBOX-COST suite.
 */
static std::string_view suite_sbox_cost_get_instances_by_year(const int year) {
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
static coco_problem_ptr_t coco_get_bbob_boxed_problem(const size_t function, const size_t dimension,
                                                      const size_t instance) {
  coco_problem_ptr_t problem =
      bbob_make_problem(function, dimension, instance, coco_bbob_variant_t::boxed, "bbob-boxed_f%03zu_i%02zu_d%02zu",
                        "bbob-boxed suite problem f%zu instance %zu in %zuD");
  return transform_vars_enforce_box(std::move(problem));
}

/**
 * @brief The bbob-boxed suite.
 */
class suite_sbox_cost_t final : public coco_suite_s {
public:
  using coco_suite_s::coco_suite_s;

  std::string_view get_instances_by_year(const int year) const override {
    return suite_sbox_cost_get_instances_by_year(year);
  }

  /** @brief Returns the problem with the given function, dimension and instance. */
  coco_problem_ptr_t get_problem(const size_t function, const size_t dimension, const size_t instance) override {
    coco_problem_ptr_t problem;

    problem = coco_get_bbob_boxed_problem(function, dimension, instance);

    return problem;
  }
};

/**
 * @brief Sets the dimensions and default instances for the bbob-boxed suite.
 */
static coco_suite_t *suite_sbox_cost_initialize() {

  coco_suite_t *suite;
  const size_t dimensions[] = {2, 3, 5, 10, 20, 40};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  suite = new suite_sbox_cost_t("bbob-boxed", 24, num_dimensions, dimensions, "year: 2023", 1);

  return suite;
}
