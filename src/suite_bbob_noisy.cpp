/**
 * @file suite_bbob_noisy.cpp
 * @brief Implementation of the bbob-noisy suite containing 30 noisy single-objective functions in 6 dimensions.
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
 * @brief Creates and returns a bbob-noisy problem without needing the actual bbob-noisy suite.
 *
 * The functions come in groups of three, which add Gaussian, uniform and Cauchy noise to the same noise-free bbob
 * problem: f101-f106 with moderate noise, f107-f130 with severe noise.
 */
static coco_problem_ptr_t coco_get_bbob_noisy_problem(const size_t function, const size_t dimension,
                                                      const size_t instance) {
  if (function < 101 || function > 130) {
    coco_error("coco_get_bbob_noisy_problem(): cannot retrieve problem f%zu instance %zu in %zuD", function, instance,
               dimension);
    return nullptr; /* Never reached */
  }

  const char *problem_id_template = "bbob_noisy_f%zu_i%02zu_d%02zu";
  const char *problem_name_template = "BBOB-NOISY suite problem f%zu instance %zu in %zuD";

  /* The seed of the bbob function the noise-free problem is based on */
  const auto rseed = [instance](const size_t bbob_function) {
    return static_cast<long>(bbob_function + 10000 * instance);
  };

  coco_problem_ptr_t problem;
  switch ((function - 101) / 3) {
  case 0: /* f101-f103 */
  case 2: /* f107-f109 */
    problem = std::make_unique<bbob_sphere_t>(function, dimension, instance, rseed(1), coco_bbob_variant_t::bbob,
                                              problem_id_template, problem_name_template);
    break;
  case 1: /* f104-f106 */
  case 3: /* f110-f112 */
    problem = std::make_unique<bbob_rosenbrock_t>(function, dimension, instance, rseed(8), coco_bbob_variant_t::bbob,
                                                  problem_id_template, problem_name_template);
    break;
  case 4: /* f113-f115 */
    problem =
        std::make_unique<bbob_step_ellipsoid_t>(function, dimension, instance, rseed(7), coco_bbob_variant_t::bbob, 0.0,
                                                problem_id_template, problem_name_template);
    break;
  case 5: /* f116-f118 */
    problem =
        std::make_unique<bbob_ellipsoid_rotated_t>(function, dimension, instance, rseed(10), coco_bbob_variant_t::bbob,
                                                   1.0e4, problem_id_template, problem_name_template);
    break;
  case 6: /* f119-f121 */
    problem =
        std::make_unique<bbob_different_powers_t>(function, dimension, instance, rseed(14), coco_bbob_variant_t::bbob,
                                                  problem_id_template, problem_name_template);
    break;
  case 7: /* f122-f124 */
    problem = std::make_unique<bbob_schaffers_t>(function, dimension, instance, rseed(17), coco_bbob_variant_t::bbob,
                                                 10.0, 0.0, problem_id_template, problem_name_template);
    break;
  case 8: /* f125-f127 */
    problem = std::make_unique<bbob_griewank_rosenbrock_t>(function, dimension, instance, rseed(19), 1.0,
                                                           problem_id_template, problem_name_template);
    break;
  default: /* f128-f130 */
    problem =
        std::make_unique<bbob_gallagher_t>(function, dimension, instance, rseed(21), coco_bbob_variant_t::bbob,
                                           static_cast<size_t>(101), 0.0, problem_id_template, problem_name_template);
    break;
  }

  const bool moderate = function <= 106;
  const double level = moderate ? 0.01 : 1.0;
  switch ((function - 101) % 3) {
  case 0:
    return transform_obj_gaussian_noise(std::move(problem), level);
  case 1:
    return transform_obj_uniform_noise(std::move(problem), level * (0.49 + 1.0 / static_cast<double>(dimension)),
                                       level);
  default:
    return transform_obj_cauchy_noise(std::move(problem), level, moderate ? 0.05 : 0.2);
  }
}

/**
 * @brief The bbob-noisy suite.
 */
class suite_bbob_noisy_t final : public coco_suite_s {
public:
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
static coco_suite_t *suite_bbob_noisy_initialize() {

  coco_suite_t *suite;
  const size_t dimensions[] = {2, 3, 5, 10, 20, 40};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  /* IMPORTANT: Make sure to change the default instance for every new workshop! */
  suite = new suite_bbob_noisy_t("bbob-noisy", 30, num_dimensions, dimensions, "year:2009", 1);

  return suite;
}
