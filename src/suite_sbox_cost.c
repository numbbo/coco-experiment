/** @file suite_sbox_cost.c
 *
 * @brief Implementation of the SBOX-COST suite, containing the 24 noiseless
 * single-objective functions in 6 dimensions from the original BBOB suite.
 *
 * FIXME: Drop f9, f19, f20 and f24 which have optima that are always close to
 *   the origin?
 */

#include "coco.h"

#include "f_attractive_sector.c"
#include "f_bent_cigar.c"
#include "f_bueche_rastrigin.c"
#include "f_different_powers.c"
#include "f_discus.c"
#include "f_ellipsoid.c"
#include "f_gallagher.c"
#include "f_griewank_rosenbrock.c"
#include "f_katsuura.c"
#include "f_linear_slope.c"
#include "f_lunacek_bi_rastrigin.c"
#include "f_rastrigin.c"
#include "f_rosenbrock.c"
#include "f_schaffers.c"
#include "f_schwefel.c"
#include "f_sharp_ridge.c"
#include "f_sphere.c"
#include "f_step_ellipsoid.c"
#include "f_weierstrass.c"

#include "transform_vars_enforce_box.c"

static coco_suite_t *coco_suite_allocate(const char *suite_name, const size_t number_of_functions,
                                         const size_t number_of_dimensions, const size_t *dimensions,
                                         const char *default_instances, const int known_optima);

/**
 * @brief Sets the dimensions and default instances for the SBOX-COST suite.
 */
static coco_suite_t *suite_sbox_cost_initialize(void) {

  coco_suite_t *suite;
  const size_t dimensions[] = {2, 3, 5, 10, 20, 40};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  suite = coco_suite_allocate("sbox-cost", 24, num_dimensions, dimensions, "year: 2023", 1);

  return suite;
}

/**
 * @brief Sets the instances associated with years for the SBOX-COST suite.
 */
static const char *suite_sbox_cost_get_instances_by_year(const int year) {
  if (year == 2023) {
    return "1-5,101-110";
  } else {
    coco_error(
        "suite_sbox_cost_get_instances_by_year(): year %d not defined for "
        "suite_sbox_cost",
        year);
    return NULL;
  }
}

/**
 * @brief Creates and returns a SBOX-COST problem without needing the actual
 * SBOX-COST suite.
 */
static coco_problem_t *coco_get_sbox_cost_problem(const size_t function,
                                                  const size_t dimension,
                                                  const size_t instance) {
  coco_problem_t *problem = NULL;

  const char *problem_id_template = "sbox-cost_f%03lu_i%02lu_d%02lu";
  const char *problem_name_template =
      "SBOX-COST suite problem f%lu instance %lu in %luD";

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
    f_step_ellipsoid_args_t args;
    args.penalty_scale = 1.0;
    problem = f_step_ellipsoid_bbob_problem_allocate(function, dimension, instance, rseed, &args, problem_id_template,
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
    f_ellipsoid_args_t args;
    args.conditioning = 1.0e6;
    problem = f_ellipsoid_rotated_bbob_problem_allocate(function, dimension, instance, rseed, &args,
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
    f_schaffers_args_t args;
    args.conditioning = 10.0;
    args.penalty_scale = 10.0;
    problem = f_schaffers_bbob_problem_allocate(function, dimension, instance, rseed, &args, problem_id_template,
                                                problem_name_template);
  } else if (function == 18) {
    f_schaffers_args_t args;
    args.conditioning = 1000.0;
    args.penalty_scale = 10.0;
    problem = f_schaffers_bbob_problem_allocate(function, dimension, instance, rseed_17, &args, problem_id_template,
                                                problem_name_template);
  } else if (function == 19) {
    f_griewank_rosenbrock_args_t args;
    args.facftrue = 10.0;
    problem = f_griewank_rosenbrock_bbob_problem_allocate(function, dimension, instance, rseed, &args,
                                                        problem_id_template, problem_name_template);
  } else if (function == 20) {
    problem = f_schwefel_bbob_problem_allocate(function, dimension, instance,
                                               rseed, problem_id_template,
                                               problem_name_template);
  } else if (function == 21) {
    f_gallagher_args_t args;
    args.number_of_peaks = (size_t)101;
    args.penalty_scale = 1.0;
    problem = f_gallagher_bbob_problem_allocate(function, dimension, instance, rseed, &args, problem_id_template,
                                                problem_name_template);
  } else if (function == 22) {
    f_gallagher_args_t args;
    args.number_of_peaks = (size_t)21;
    args.penalty_scale = 1.0;
    problem = f_gallagher_bbob_problem_allocate(function, dimension, instance, rseed, &args, problem_id_template,
                                                problem_name_template);
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
        "coco_get_sbox_cost_problem(): cannot retrieve problem f%lu instance "
        "%lu in %luD",
        (unsigned long)function, (unsigned long)instance,
        (unsigned long)dimension);
    return NULL; /* Never reached */
  }

  problem = transform_vars_enforce_box(problem);
  return problem;
}

/**
 * @brief Returns the problem from the SBOX-COST suite that corresponds to the
 * given parameters.
 *
 * @param suite The COCO suite.
 * @param function_idx Index of the function (starting from 0).
 * @param dimension_idx Index of the dimension (starting from 0).
 * @param instance_idx Index of the instance (starting from 0).
 * @return The problem that corresponds to the given parameters.
 */
static coco_problem_t *suite_sbox_cost_get_problem(coco_suite_t *suite,
                                                   const size_t function_idx,
                                                   const size_t dimension_idx,
                                                   const size_t instance_idx) {

  coco_problem_t *problem = NULL;

  const size_t function = suite->functions[function_idx];
  const size_t dimension = suite->dimensions[dimension_idx];
  const size_t instance = suite->instances[instance_idx];

  problem = coco_get_sbox_cost_problem(function, dimension, instance);

  problem->suite_dep_function = function;
  problem->suite_dep_instance = instance;
  problem->suite_dep_index = coco_suite_encode_problem_index(
      suite, function_idx, dimension_idx, instance_idx);

  return problem;
}
