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

#include "f_attractive_sector.cpp"
#include "f_bent_cigar.cpp"
#include "f_bueche_rastrigin.cpp"
#include "f_different_powers.cpp"
#include "f_discus.cpp"
#include "f_ellipsoid.cpp"
#include "f_gallagher.cpp"
#include "f_griewank_rosenbrock.cpp"
#include "f_katsuura.cpp"
#include "f_linear_slope.cpp"
#include "f_lunacek_bi_rastrigin.cpp"
#include "f_rastrigin.cpp"
#include "f_rosenbrock.cpp"
#include "f_schaffers.cpp"
#include "f_schwefel.cpp"
#include "f_sharp_ridge.cpp"
#include "f_sphere.cpp"
#include "f_step_ellipsoid.cpp"
#include "f_weierstrass.cpp"

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
  coco_problem_ptr_t problem;

  const char *problem_id_template = "bbob-boxed_f%03lu_i%02lu_d%02lu";
  const char *problem_name_template =
      "bbob-boxed suite problem f%lu instance %lu in %luD";

  const long rseed = (long)(function + 10000 * instance);
  const long rseed_3 = (long)(3 + 10000 * instance);
  const long rseed_17 = (long)(17 + 10000 * instance);

  if (function == 1) {
    problem = f_sphere_bbob_problem_allocate(function, dimension, instance,
                                             rseed, problem_id_template,
                                             problem_name_template);
  } else if (function == 2) {
    problem = f_ellipsoid_bbob_problem_allocate(function, dimension, instance,
                                                rseed, problem_id_template,
                                                problem_name_template);
  } else if (function == 3) {
    problem = f_rastrigin_bbob_problem_allocate(function, dimension, instance,
                                                rseed, problem_id_template,
                                                problem_name_template);
  } else if (function == 4) {
    problem = f_bueche_rastrigin_bbob_problem_allocate(
        function, dimension, instance, rseed_3, problem_id_template,
        problem_name_template);
  } else if (function == 5) {
    problem = f_linear_slope_bbob_problem_allocate(
        function, dimension, instance, rseed, problem_id_template,
        problem_name_template);
  } else if (function == 6) {
    problem = f_attractive_sector_bbob_problem_allocate(
        function, dimension, instance, rseed, problem_id_template,
        problem_name_template);
  } else if (function == 7) {
    problem = f_step_ellipsoid_bbob_problem_allocate(function, dimension, instance, rseed, 1.0, problem_id_template,
                                                     problem_name_template);
  } else if (function == 8) {
    problem = f_rosenbrock_bbob_problem_allocate(function, dimension, instance,
                                                 rseed, problem_id_template,
                                                 problem_name_template);
  } else if (function == 9) {
    problem = f_rosenbrock_rotated_bbob_problem_allocate(
        function, dimension, instance, rseed, problem_id_template,
        problem_name_template);
  } else if (function == 10) {
    problem = f_ellipsoid_rotated_bbob_problem_allocate(function, dimension, instance, rseed, 1.0e6,
                                                        problem_id_template, problem_name_template);
  } else if (function == 11) {
    problem = f_discus_bbob_problem_allocate(function, dimension, instance,
                                             rseed, problem_id_template,
                                             problem_name_template);
  } else if (function == 12) {
    problem = f_bent_cigar_bbob_problem_allocate(function, dimension, instance,
                                                 rseed, problem_id_template,
                                                 problem_name_template);
  } else if (function == 13) {
    problem = f_sharp_ridge_bbob_problem_allocate(function, dimension, instance,
                                                  rseed, problem_id_template,
                                                  problem_name_template);
  } else if (function == 14) {
    problem = f_different_powers_bbob_problem_allocate(
        function, dimension, instance, rseed, problem_id_template,
        problem_name_template);
  } else if (function == 15) {
    problem = f_rastrigin_rotated_bbob_problem_allocate(
        function, dimension, instance, rseed, problem_id_template,
        problem_name_template);
  } else if (function == 16) {
    problem = f_weierstrass_bbob_problem_allocate(function, dimension, instance,
                                                  rseed, problem_id_template,
                                                  problem_name_template);
  } else if (function == 17) {
    problem = f_schaffers_bbob_problem_allocate(function, dimension, instance, rseed, 10.0, 10.0, problem_id_template,
                                                problem_name_template);
  } else if (function == 18) {
    problem = f_schaffers_bbob_problem_allocate(function, dimension, instance, rseed_17, 1000.0, 10.0,
                                                problem_id_template, problem_name_template);
  } else if (function == 19) {
    problem = f_griewank_rosenbrock_bbob_problem_allocate(function, dimension, instance, rseed, 10.0,
                                                          problem_id_template, problem_name_template);
  } else if (function == 20) {
    problem = f_schwefel_bbob_problem_allocate(function, dimension, instance,
                                               rseed, problem_id_template,
                                               problem_name_template);
  } else if (function == 21) {
    problem = f_gallagher_bbob_problem_allocate(function, dimension, instance, rseed, (size_t)101, 1.0,
                                                problem_id_template, problem_name_template);
  } else if (function == 22) {
    problem = f_gallagher_bbob_problem_allocate(function, dimension, instance, rseed, (size_t)21, 1.0,
                                                problem_id_template, problem_name_template);
  } else if (function == 23) {
    problem = f_katsuura_bbob_problem_allocate(function, dimension, instance,
                                               rseed, problem_id_template,
                                               problem_name_template);
  } else if (function == 24) {
    problem = f_lunacek_bi_rastrigin_bbob_problem_allocate(
        function, dimension, instance, rseed, problem_id_template,
        problem_name_template);
  } else {
    coco_error(
        "coco_get_bbob_boxed_problem(): cannot retrieve problem f%lu instance "
        "%lu in %luD",
        (unsigned long)function, (unsigned long)instance,
        (unsigned long)dimension);
    return nullptr; /* Never reached */
  }

  problem = transform_vars_enforce_box(std::move(problem));
  return problem;
}

/**
 * @brief The bbob-boxed suite.
 */
struct suite_sbox_cost_t final : public coco_suite_s {
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
static coco_suite_t *suite_sbox_cost_initialize(void) {

  coco_suite_t *suite;
  const size_t dimensions[] = {2, 3, 5, 10, 20, 40};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  suite = new suite_sbox_cost_t("bbob-boxed", 24, num_dimensions, dimensions, "year: 2023", 1);

  return suite;
}
