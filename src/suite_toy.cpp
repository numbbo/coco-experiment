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

namespace coco {

/**
 * @brief The toy suite.
 */
class SuiteToy final : public Suite {
public:
  using Suite::Suite;

  /** @brief Returns the problem with the given function, dimension and instance. */
  ProblemPtr GetProblem(const size_t function, const size_t dimension,
                        const size_t /* instance: the toy functions have none */) override {
    ProblemPtr problem;

    if (function == 1) {
      problem = std::make_unique<Sphere>(dimension);
    } else if (function == 2) {
      problem = std::make_unique<Ellipsoid>(dimension, 1.0e6);
    } else if (function == 3) {
      problem = std::make_unique<Rastrigin>(dimension);
    } else if (function == 4) {
      problem = std::make_unique<BuecheRastrigin>(dimension);
    } else if (function == 5) {
      double xopt[40] = {5.0};
      problem = std::make_unique<LinearSlope>(dimension, std::span(xopt).first(dimension));
    } else if (function == 6) {
      problem = std::make_unique<Rosenbrock>(dimension);
    } else {
      coco_error("suite_toy_t::GetProblem(): function %zu does not exist in this suite", function);
      return nullptr; /* Never reached */
    }

    return problem;
  }
};

/**
 * @brief Sets the dimensions and default instances for the toy suite.
 */
static Suite *SuiteToyInitialize() {

  Suite *suite;
  const size_t dimensions[] = {2, 3, 5, 10, 20};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  suite = new SuiteToy("toy", 6, num_dimensions, dimensions, "instances: 1", 1);

  return suite;
}

} // namespace coco
