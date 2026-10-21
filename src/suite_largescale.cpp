/**
 * @file suite_largescale.cpp
 * @brief Implementation of the bbob large-scale suite containing 24 functions in 6 large dimensions.
 */

#include "coco.h"

#include "f_attractive_sector.cpp"
#include "f_bent_cigar_generalized.cpp"
#include "f_bueche_rastrigin.cpp"
#include "f_different_powers.cpp"
#include "f_discus_generalized.cpp"
#include "f_ellipsoid.cpp"
#include "f_griewank_rosenbrock.cpp"
#include "f_linear_slope.cpp"
#include "f_rastrigin.cpp"
#include "f_schaffers.cpp"
#include "f_schwefel_generalized.cpp"
#include "f_sharp_ridge_generalized.cpp"
#include "f_sphere.cpp"
#include "f_step_ellipsoid.cpp"
#include "f_weierstrass.cpp"

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
  coco_problem_ptr_t problem;

  const char *problem_id_template = "bbob_f%03lu_i%02lu_d%04lu";
  const char *problem_name_template = "BBOB large-scale suite problem f%lu instance %lu in %luD";

  const long rseed = (long)(function + 10000 * instance);
  const long rseed_3 = (long)(3 + 10000 * instance);
  const long rseed_17 = (long)(17 + 10000 * instance);

  if (function == 1) {
    problem = f_sphere_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::largescale,
                                             problem_id_template, problem_name_template);
  } else if (function == 2) {
    problem = f_ellipsoid_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::largescale,
                                                problem_id_template, problem_name_template);
  } else if (function == 3) {
    problem = f_rastrigin_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::largescale,
                                                problem_id_template, problem_name_template);
  } else if (function == 4) {
    problem = f_bueche_rastrigin_bbob_problem_allocate(function, dimension, instance, rseed_3,
                                                       coco_bbob_variant_t::largescale, problem_id_template,
                                                       problem_name_template);
  } else if (function == 5) {
    problem =
        f_linear_slope_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::largescale,
                                             problem_id_template, problem_name_template);
  } else if (function == 6) {
    problem = f_attractive_sector_permblockdiag_bbob_problem_allocate(function, dimension, instance, rseed,
                                                                      coco_bbob_variant_t::largescale,
                                                                      problem_id_template, problem_name_template);
  } else if (function == 7) {
    problem = f_step_ellipsoid_permblockdiag_bbob_problem_allocate(function, dimension, instance, rseed,
                                                                   coco_bbob_variant_t::largescale, problem_id_template,
                                                                   problem_name_template);
  } else if (function == 8) {
    problem = f_rosenbrock_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::largescale,
                                                 problem_id_template, problem_name_template);
  } else if (function == 9) {
    problem = f_rosenbrock_permblockdiag_bbob_problem_allocate(function, dimension, instance, rseed,
                                                               coco_bbob_variant_t::largescale, problem_id_template,
                                                               problem_name_template);
  } else if (function == 10) {
    problem = f_ellipsoid_permblockdiag_bbob_problem_allocate(function, dimension, instance, rseed,
                                                              coco_bbob_variant_t::largescale, problem_id_template,
                                                              problem_name_template);
  } else if (function == 11) {
    problem = f_discus_generalized_permblockdiag_bbob_problem_allocate(function, dimension, instance, rseed,
                                                                       coco_bbob_variant_t::largescale,
                                                                       problem_id_template, problem_name_template);
  } else if (function == 12) {
    problem = f_bent_cigar_generalized_permblockdiag_bbob_problem_allocate(function, dimension, instance, rseed,
                                                                           coco_bbob_variant_t::largescale,
                                                                           problem_id_template, problem_name_template);
  } else if (function == 13) {
    problem = f_sharp_ridge_generalized_permblockdiag_bbob_problem_allocate(function, dimension, instance, rseed,
                                                                            coco_bbob_variant_t::largescale,
                                                                            problem_id_template, problem_name_template);
  } else if (function == 14) {
    problem = f_different_powers_permblockdiag_bbob_problem_allocate(function, dimension, instance, rseed,
                                                                     coco_bbob_variant_t::largescale,
                                                                     problem_id_template, problem_name_template);
  } else if (function == 15) {
    problem = f_rastrigin_permblockdiag_bbob_problem_allocate(function, dimension, instance, rseed,
                                                              coco_bbob_variant_t::largescale, problem_id_template,
                                                              problem_name_template);
  } else if (function == 16) {
    problem = f_weierstrass_permblockdiag_bbob_problem_allocate(function, dimension, instance, rseed,
                                                                coco_bbob_variant_t::largescale, problem_id_template,
                                                                problem_name_template);
  } else if (function == 17) {
    problem = f_schaffers_permblockdiag_bbob_problem_allocate(function, dimension, instance, rseed,
                                                              coco_bbob_variant_t::largescale, 10, problem_id_template,
                                                              problem_name_template);
  } else if (function == 18) {
    problem = f_schaffers_permblockdiag_bbob_problem_allocate(function, dimension, instance, rseed_17,
                                                              coco_bbob_variant_t::largescale, 1000,
                                                              problem_id_template, problem_name_template);
  } else if (function == 19) {
    problem = f_griewank_rosenbrock_permblockdiag_bbob_bbob_problem_allocate(
        function, dimension, instance, rseed, 10.0, problem_id_template, problem_name_template);
  } else if (function == 20) {
    problem = f_schwefel_generalized_bbob_problem_allocate(function, dimension, instance, rseed, problem_id_template,
                                                           problem_name_template);
  } else if (function == 21) {
    problem = f_gallagher_permblockdiag_bbob_problem_allocate(function, dimension, instance, rseed, 101,
                                                              problem_id_template, problem_name_template);
  } else if (function == 22) {
    problem = f_gallagher_permblockdiag_bbob_problem_allocate(function, dimension, instance, rseed, 21,
                                                              problem_id_template, problem_name_template);
  } else if (function == 23) {
    problem = f_katsuura_permblockdiag_bbob_problem_allocate(function, dimension, instance, rseed,
                                                             coco_bbob_variant_t::largescale, problem_id_template,
                                                             problem_name_template);
  } else if (function == 24) {
    problem = f_lunacek_bi_rastrigin_permblockdiag_bbob_problem_allocate(function, dimension, instance, rseed,
                                                                         problem_id_template, problem_name_template);
  } else {
    coco_error("coco_get_largescale_problem(): cannot retrieve problem f%lu instance %lu in %luD",
               (unsigned long)function, (unsigned long)instance, (unsigned long)dimension);
    return nullptr; /* Never reached */
  }

  if (problem == nullptr) {
    coco_error("coco_get_largescale_problem(): function f%lu not yet implemented  ", function);
  }
  return problem;
}

/**
 * @brief The bbob-largescale suite.
 */
struct suite_largescale_t final : public coco_suite_s {
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
static coco_suite_t *suite_largescale_initialize(void) {

  coco_suite_t *suite;
  const size_t dimensions[] = {20, 40, 80, 160, 320, 640};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  suite = new suite_largescale_t("bbob-largescale", 24, num_dimensions, dimensions, "instances: 1-15", 1);
  return suite;
}
