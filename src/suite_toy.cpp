/**
 * @file suite_toy.cpp
 * @brief Implementation of a toy suite containing 6 noiseless "basic" single-objective functions in 5
 * dimensions.
 */

#include "coco_problem.cpp"
#include "f_bueche_rastrigin.cpp"
#include "f_ellipsoid.cpp"
#include "f_linear_slope.cpp"
#include "f_rastrigin.cpp"
#include "f_rosenbrock.cpp"
#include "f_sphere.cpp"

/**
 * @brief The toy suite.
 */
class suite_toy_t final : public coco_suite_s {
public:
  using coco_suite_s::coco_suite_s;

  /** @brief Returns the problem with the given function, dimension and instance. */
  coco_problem_ptr_t get_problem(const size_t function, const size_t dimension,
                                 const size_t /* instance: the toy functions have none */) override {
    coco_problem_ptr_t problem;

    if (function == 1) {
      problem = std::make_unique<f_sphere_t>(dimension);
    } else if (function == 2) {
      problem = std::make_unique<f_ellipsoid_t>(dimension, 1.0e6);
    } else if (function == 3) {
      problem = std::make_unique<f_rastrigin_t>(dimension);
    } else if (function == 4) {
      problem = std::make_unique<f_bueche_rastrigin_t>(dimension);
    } else if (function == 5) {
      double xopt[40] = {5.0};
      problem = std::make_unique<f_linear_slope_t>(dimension, std::span(xopt).first(dimension));
    } else if (function == 6) {
      problem = std::make_unique<f_rosenbrock_t>(dimension);
    } else {
      coco_error("suite_toy_t::get_problem(): function %zu does not exist in this suite", function);
      return nullptr; /* Never reached */
    }

    return problem;
  }
};

/**
 * @brief Sets the dimensions and default instances for the toy suite.
 */
static coco_suite_t *suite_toy_initialize() {

  coco_suite_t *suite;
  const size_t dimensions[] = {2, 3, 5, 10, 20};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  suite = new suite_toy_t("toy", 6, num_dimensions, dimensions, "instances: 1", 1);

  return suite;
}
