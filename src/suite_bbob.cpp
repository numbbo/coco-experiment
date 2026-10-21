/**
 * @file suite_bbob.cpp
 * @brief Implementation of the bbob suite containing 24 noiseless single-objective functions in 6
 * dimensions.
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

/**
 * @brief Sets the instances associated with years for the bbob suite.
 */
static std::string_view suite_bbob_get_instances_by_year(const int year) {

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
    coco_error("suite_bbob_get_instances_by_year(): year %d not defined for suite_bbob", year);
    return {};
  }
}

/**
 * @brief Creates and returns a BBOB problem without needing the actual bbob suite.
 *
 * Useful for other suites as well (see for example suite_biobj.cpp).
 */
static coco_problem_ptr_t coco_make_bbob_problem(const size_t function, const size_t dimension, const size_t instance) {
  coco_problem_ptr_t problem;

  const char *problem_id_template = "bbob_f%03lu_i%02lu_d%02lu";
  const char *problem_name_template = "BBOB suite problem f%lu instance %lu in %luD";

  const long rseed = (long)(function + 10000 * instance);
  const long rseed_3 = (long)(3 + 10000 * instance);
  const long rseed_17 = (long)(17 + 10000 * instance);

  if (function == 1) {
    problem = f_sphere_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                             problem_id_template, problem_name_template);
  } else if (function == 2) {
    problem = f_ellipsoid_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                problem_id_template, problem_name_template);
  } else if (function == 3) {
    problem = f_rastrigin_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                problem_id_template, problem_name_template);
  } else if (function == 4) {
    problem = f_bueche_rastrigin_bbob_problem_allocate(
        function, dimension, instance, rseed_3, coco_bbob_variant_t::bbob, problem_id_template, problem_name_template);
  } else if (function == 5) {
    problem = f_linear_slope_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                   problem_id_template, problem_name_template);
  } else if (function == 6) {
    problem = f_attractive_sector_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                        problem_id_template, problem_name_template);
  } else if (function == 7) {
    problem = f_step_ellipsoid_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                     1.0, problem_id_template, problem_name_template);
  } else if (function == 8) {
    problem = f_rosenbrock_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                 problem_id_template, problem_name_template);
  } else if (function == 9) {
    problem = f_rosenbrock_rotated_bbob_problem_allocate(function, dimension, instance, rseed, problem_id_template,
                                                         problem_name_template);
  } else if (function == 10) {
    problem = f_ellipsoid_rotated_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                        1.0e6, problem_id_template, problem_name_template);
  } else if (function == 11) {
    problem = f_discus_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                             problem_id_template, problem_name_template);
  } else if (function == 12) {
    problem = f_bent_cigar_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                 problem_id_template, problem_name_template);
  } else if (function == 13) {
    problem = f_sharp_ridge_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                  problem_id_template, problem_name_template);
  } else if (function == 14) {
    problem = f_different_powers_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                       problem_id_template, problem_name_template);
  } else if (function == 15) {
    problem = f_rastrigin_rotated_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                        problem_id_template, problem_name_template);
  } else if (function == 16) {
    problem = f_weierstrass_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                  problem_id_template, problem_name_template);
  } else if (function == 17) {
    problem = f_schaffers_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob, 10.0,
                                                10.0, problem_id_template, problem_name_template);
  } else if (function == 18) {
    problem = f_schaffers_bbob_problem_allocate(function, dimension, instance, rseed_17, coco_bbob_variant_t::bbob,
                                                1000.0, 10.0, problem_id_template, problem_name_template);
  } else if (function == 19) {
    problem = f_griewank_rosenbrock_bbob_problem_allocate(function, dimension, instance, rseed, 10.0,
                                                          problem_id_template, problem_name_template);
  } else if (function == 20) {
    problem = f_schwefel_bbob_problem_allocate(function, dimension, instance, rseed, problem_id_template,
                                               problem_name_template);
  } else if (function == 21) {
    problem = f_gallagher_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                (size_t)101, 1.0, problem_id_template, problem_name_template);
  } else if (function == 22) {
    problem = f_gallagher_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                (size_t)21, 1.0, problem_id_template, problem_name_template);
  } else if (function == 23) {
    problem = f_katsuura_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                               problem_id_template, problem_name_template);
  } else if (function == 24) {
    problem = f_lunacek_bi_rastrigin_bbob_problem_allocate(
        function, dimension, instance, rseed, coco_bbob_variant_t::bbob, problem_id_template, problem_name_template);
  } else {
    coco_error("coco_get_bbob_problem(): cannot retrieve problem f%lu instance %lu in %luD", (unsigned long)function,
               (unsigned long)instance, (unsigned long)dimension);
    return nullptr; /* Never reached */
  }

  return problem;
}

/**
 * @brief Creates and returns a BBOB problem, which the caller must free with coco_problem_free().
 *
 * Used by the language bindings.
 */
coco_problem_t *coco_get_bbob_problem(const size_t function, const size_t dimension, const size_t instance) {
  return coco_make_bbob_problem(function, dimension, instance).release();
}

/**
 * @brief The bbob suite.
 */
struct suite_bbob_t final : public coco_suite_s {
  using coco_suite_s::coco_suite_s;

  std::string_view get_instances_by_year(const int year) const override {
    return suite_bbob_get_instances_by_year(year);
  }

  /** @brief Returns the problem with the given function, dimension and instance. */
  coco_problem_ptr_t get_problem(const size_t function, const size_t dimension, const size_t instance) override {
    coco_problem_ptr_t problem;

    problem = coco_make_bbob_problem(function, dimension, instance);

    return problem;
  }
};

/**
 * @brief Sets the dimensions and default instances for the bbob suite.
 */
static coco_suite_t *suite_bbob_initialize(void) {

  coco_suite_t *suite;
  const size_t dimensions[] = {2, 3, 5, 10, 20, 40};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  /* IMPORTANT: Make sure to change the default instance for every new workshop! */
  suite = new suite_bbob_t("bbob", 24, num_dimensions, dimensions, "year: 2018", 1);

  return suite;
}
