/**
 * @file suite_bbob_noisy.cpp
 * @brief Implementation of the bbob suite containing 24 noiseless single-objective functions in 6
 * dimensions.
 */
#include "transform_obj_gaussian_noise.cpp"
#include "transform_obj_uniform_noise.cpp"
#include "transform_obj_cauchy_noise.cpp"

#include "f_different_powers.cpp"
#include "f_ellipsoid.cpp"
#include "f_gallagher.cpp"
#include "f_griewank_rosenbrock.cpp"
#include "f_rosenbrock.cpp"
#include "f_schaffers.cpp"
#include "f_sphere.cpp"
#include "f_step_ellipsoid.cpp"

/**
 * @brief Sets the instances associated with years for the bbob suite.
 */
static std::string_view suite_bbob_noisy_get_instances_by_year(const int year) {

  if (year <= 2009) {
    return "1-15";
  } else {
    coco_error("suite_bbob_noisy_get_instances_by_year(): year %d not defined for suite bbob-noisy", year);
    return {};
  }
}

/**
 * @brief Creates and returns a BBOB problem without needing the actual bbob suite.
 *
 * Useful for other suites as well (see for example suite_biobj.cpp).
 */
static coco_problem_ptr_t coco_get_bbob_noisy_problem(const size_t function, const size_t dimension,
                                                      const size_t instance) {
  coco_problem_ptr_t problem;
  coco_problem_ptr_t inner_problem;

  const char *problem_id_template = "bbob_noisy_f%lu_i%02lu_d%02lu";
  const char *problem_name_template = "BBOB-NOISY suite problem f%lu instance %lu in %luD";

  long rseed;
  const long rseed_1 = (long)(1 + 10000 * instance);
  const long rseed_8 = (long)(8 + 10000 * instance);
  const long rseed_7 = (long)(7 + 10000 * instance);
  const long rseed_10 = (long)(10 + 10000 * instance);
  const long rseed_14 = (long)(14 + 10000 * instance);
  const long rseed_17 = (long)(17 + 10000 * instance);
  const long rseed_19 = (long)(19 + 10000 * instance);
  const long rseed_21 = (long)(21 + 10000 * instance);

  if (function == 101) {
    rseed = rseed_1;
    inner_problem = f_sphere_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                   problem_id_template, problem_name_template);
    double beta = 0.01;
    problem = transform_obj_gaussian_noise(std::move(inner_problem), beta);
  } else if (function == 102) {
    rseed = rseed_1;
    inner_problem = f_sphere_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                   problem_id_template, problem_name_template);
    double alpha = 0.01 * (0.49 + 1 / (double)dimension);
    double beta = 0.01;
    problem = transform_obj_uniform_noise(std::move(inner_problem), alpha, beta);
  } else if (function == 103) {
    rseed = rseed_1;
    inner_problem = f_sphere_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                   problem_id_template, problem_name_template);
    double alpha = 0.01;
    double p = 0.05;
    problem = transform_obj_cauchy_noise(std::move(inner_problem), alpha, p);
  } else if (function == 104) {
    rseed = rseed_8;
    inner_problem = f_rosenbrock_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                       problem_id_template, problem_name_template);
    double beta = 0.01;
    problem = transform_obj_gaussian_noise(std::move(inner_problem), beta);
  } else if (function == 105) {
    rseed = rseed_8;
    inner_problem = f_rosenbrock_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                       problem_id_template, problem_name_template);
    double alpha = 0.01 * (0.49 + 1 / (double)dimension);
    double beta = 0.01;
    problem = transform_obj_uniform_noise(std::move(inner_problem), alpha, beta);
  } else if (function == 106) {
    rseed = rseed_8;
    inner_problem = f_rosenbrock_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                       problem_id_template, problem_name_template);
    double alpha = 0.01;
    double p = 0.05;
    problem = transform_obj_cauchy_noise(std::move(inner_problem), alpha, p);
  } else if (function == 107) {
    rseed = rseed_1;
    inner_problem = f_sphere_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                   problem_id_template, problem_name_template);
    double beta = 1.0;
    problem = transform_obj_gaussian_noise(std::move(inner_problem), beta);
  } else if (function == 108) {
    rseed = rseed_1;
    inner_problem = f_sphere_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                   problem_id_template, problem_name_template);
    double alpha = 0.49 + 1.0 / (double)dimension;
    double beta = 1.0;
    problem = transform_obj_uniform_noise(std::move(inner_problem), alpha, beta);
  } else if (function == 109) {
    rseed = rseed_1;
    inner_problem = f_sphere_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                   problem_id_template, problem_name_template);
    double alpha = 1.0;
    double p = 0.2;
    problem = transform_obj_cauchy_noise(std::move(inner_problem), alpha, p);
  } else if (function == 110) {
    rseed = rseed_8;
    inner_problem = f_rosenbrock_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                       problem_id_template, problem_name_template);
    double beta = 1.0;
    problem = transform_obj_gaussian_noise(std::move(inner_problem), beta);
  } else if (function == 111) {
    rseed = rseed_8;
    inner_problem = f_rosenbrock_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                       problem_id_template, problem_name_template);
    double alpha = 0.49 + 1.0 / (double)dimension;
    double beta = 1.0;
    problem = transform_obj_uniform_noise(std::move(inner_problem), alpha, beta);
  } else if (function == 112) {
    rseed = rseed_8;
    inner_problem = f_rosenbrock_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                       problem_id_template, problem_name_template);
    double alpha = 1.0;
    double p = 0.2;
    problem = transform_obj_cauchy_noise(std::move(inner_problem), alpha, p);
  } else if (function == 113) {
    rseed = rseed_7;
    inner_problem =
        f_step_ellipsoid_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob, 0.0,
                                               problem_id_template, problem_name_template);
    double beta = 1.0;
    problem = transform_obj_gaussian_noise(std::move(inner_problem), beta);
  } else if (function == 114) {
    rseed = rseed_7;
    inner_problem =
        f_step_ellipsoid_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob, 0.0,
                                               problem_id_template, problem_name_template);
    double alpha = 0.49 + 1.0 / (double)dimension;
    double beta = 1.0;
    problem = transform_obj_uniform_noise(std::move(inner_problem), alpha, beta);
  } else if (function == 115) {
    rseed = rseed_7;
    inner_problem =
        f_step_ellipsoid_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob, 0.0,
                                               problem_id_template, problem_name_template);
    double alpha = 1.0;
    double p = 0.2;
    problem = transform_obj_cauchy_noise(std::move(inner_problem), alpha, p);
  } else if (function == 116) {
    rseed = rseed_10;
    inner_problem =
        f_ellipsoid_rotated_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                  1.0e4, problem_id_template, problem_name_template);
    double beta = 1.0;
    problem = transform_obj_gaussian_noise(std::move(inner_problem), beta);
  } else if (function == 117) {
    rseed = rseed_10;
    inner_problem =
        f_ellipsoid_rotated_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                  1.0e4, problem_id_template, problem_name_template);
    double alpha = 0.49 + 1.0 / (double)dimension;
    double beta = 1.0;
    problem = transform_obj_uniform_noise(std::move(inner_problem), alpha, beta);
  } else if (function == 118) {
    rseed = rseed_10;
    inner_problem =
        f_ellipsoid_rotated_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                  1.0e4, problem_id_template, problem_name_template);
    double alpha = 1.0;
    double p = 0.2;
    problem = transform_obj_cauchy_noise(std::move(inner_problem), alpha, p);
  } else if (function == 119) {
    rseed = rseed_14;
    inner_problem = f_different_powers_bbob_problem_allocate(
        function, dimension, instance, rseed, coco_bbob_variant_t::bbob, problem_id_template, problem_name_template);
    double beta = 1.0;
    problem = transform_obj_gaussian_noise(std::move(inner_problem), beta);
  } else if (function == 120) {
    rseed = rseed_14;
    inner_problem = f_different_powers_bbob_problem_allocate(
        function, dimension, instance, rseed, coco_bbob_variant_t::bbob, problem_id_template, problem_name_template);
    double alpha = 0.49 + 1.0 / (double)dimension;
    double beta = 1.0;
    problem = transform_obj_uniform_noise(std::move(inner_problem), alpha, beta);
  } else if (function == 121) {
    rseed = rseed_14;
    inner_problem = f_different_powers_bbob_problem_allocate(
        function, dimension, instance, rseed, coco_bbob_variant_t::bbob, problem_id_template, problem_name_template);
    double alpha = 1.0;
    double p = 0.2;
    problem = transform_obj_cauchy_noise(std::move(inner_problem), alpha, p);
  } else if (function == 122) {
    rseed = rseed_17;
    inner_problem = f_schaffers_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                      10.0, 0.0, problem_id_template, problem_name_template);
    double beta = 1.0;
    problem = transform_obj_gaussian_noise(std::move(inner_problem), beta);
  } else if (function == 123) {
    rseed = rseed_17;
    inner_problem = f_schaffers_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                      10.0, 0.0, problem_id_template, problem_name_template);
    double alpha = 0.49 + 1.0 / (double)dimension;
    double beta = 1.0;
    problem = transform_obj_uniform_noise(std::move(inner_problem), alpha, beta);
  } else if (function == 124) {
    rseed = rseed_17;
    inner_problem = f_schaffers_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                      10.0, 0.0, problem_id_template, problem_name_template);
    double alpha = 1.0;
    double p = 0.2;
    problem = transform_obj_cauchy_noise(std::move(inner_problem), alpha, p);
  } else if (function == 125) {
    rseed = rseed_19;
    inner_problem = f_griewank_rosenbrock_bbob_problem_allocate(function, dimension, instance, rseed, 1.0,
                                                                problem_id_template, problem_name_template);
    double beta = 1.0;
    problem = transform_obj_gaussian_noise(std::move(inner_problem), beta);
  } else if (function == 126) {
    rseed = rseed_19;
    inner_problem = f_griewank_rosenbrock_bbob_problem_allocate(function, dimension, instance, rseed, 1.0,
                                                                problem_id_template, problem_name_template);
    double alpha = 0.49 + 1.0 / (double)dimension;
    double beta = 1.0;
    problem = transform_obj_uniform_noise(std::move(inner_problem), alpha, beta);
  } else if (function == 127) {
    rseed = rseed_19;
    inner_problem = f_griewank_rosenbrock_bbob_problem_allocate(function, dimension, instance, rseed, 1.0,
                                                                problem_id_template, problem_name_template);
    double alpha = 1.0;
    double p = 0.2;
    problem = transform_obj_cauchy_noise(std::move(inner_problem), alpha, p);
  } else if (function == 128) {
    rseed = rseed_21;
    inner_problem = f_gallagher_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                      (size_t)101, 0.0, problem_id_template, problem_name_template);
    double beta = 1.0;
    problem = transform_obj_gaussian_noise(std::move(inner_problem), beta);
  } else if (function == 129) {
    rseed = rseed_21;
    inner_problem = f_gallagher_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                      (size_t)101, 0.0, problem_id_template, problem_name_template);
    double alpha = 0.49 + 1.0 / (double)dimension;
    double beta = 1.0;
    problem = transform_obj_uniform_noise(std::move(inner_problem), alpha, beta);
  } else if (function == 130) {
    rseed = rseed_21;
    inner_problem = f_gallagher_bbob_problem_allocate(function, dimension, instance, rseed, coco_bbob_variant_t::bbob,
                                                      (size_t)101, 0.0, problem_id_template, problem_name_template);
    double alpha = 1.0;
    double p = 0.2;
    problem = transform_obj_cauchy_noise(std::move(inner_problem), alpha, p);
  } else {
    coco_error("coco_get_bbob_noisy_problem(): cannot retrieve problem f%lu instance %lu in %luD",
               (unsigned long)function, (unsigned long)instance, (unsigned long)dimension);
    return nullptr; /* Never reached */
  }
  return problem;
}

/**
 * @brief The bbob-noisy suite.
 */
struct suite_bbob_noisy_t final : public coco_suite_s {
  using coco_suite_s::coco_suite_s;

  std::string_view get_instances_by_year(const int year) const override {
    return suite_bbob_noisy_get_instances_by_year(year);
  }

  /** @brief Returns the problem with the given function, dimension and instance. */
  coco_problem_ptr_t get_problem(const size_t function, const size_t dimension, const size_t instance) override {
    coco_problem_ptr_t problem;

    problem = coco_get_bbob_noisy_problem(function, dimension, instance);

    return problem;
  }
};

/**
 * @brief Sets the dimensions and default instances for the bbob suite.
 */
static coco_suite_t *suite_bbob_noisy_initialize(void) {

  coco_suite_t *suite;
  const size_t dimensions[] = {2, 3, 5, 10, 20, 40};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  /* IMPORTANT: Make sure to change the default instance for every new workshop! */
  suite = new suite_bbob_noisy_t("bbob-noisy", 30, num_dimensions, dimensions, "year:2009", 1);

  return suite;
}
