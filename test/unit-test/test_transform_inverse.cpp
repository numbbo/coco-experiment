#include "minunit.h"

#include "coco.cpp"

namespace coco {

static const size_t DIMENSION = 6;

/**
 * A base problem that remembers the decision values it was last evaluated at.
 */
struct Recorder final : public Problem {
  std::vector<double> last_x;

  Recorder() : Problem(DIMENSION, 1, 0), last_x(DIMENSION) {
    for (size_t i = 0; i < DIMENSION; ++i) {
      smallest_values_of_interest_[i] = -5.0;
      largest_values_of_interest_[i] = 5.0;
      best_parameter_[i] = 0.0;
    }
    best_value_[0] = 0.0;
    SetId("recorder_d%02lu", (unsigned long)DIMENSION);
  }

  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    last_x.assign(x.begin(), x.end());
    y[0] = 0.0;
  }
};

/** Some point in the coordinates of the base problem, with both signs and different magnitudes. */
static std::vector<double> base_point() {
  std::vector<double> base_x(DIMENSION);
  for (size_t i = 0; i < DIMENSION; ++i)
    base_x[i] = 3.0 * sin(1.0 + (double)i) * (i % 3 == 0 ? 0.1 : 1.0);
  return base_x;
}

/**
 * Checks that evaluating the problem built by make at FromBaseCoordinates(base_x) passes base_x to the base
 * problem again, up to the given relative precision (0 means exactly).
 */
template <typename Make> static int maps_back(Make make, const double precision) {
  auto owned = std::make_unique<Recorder>();
  Recorder *recorder = owned.get();
  ProblemPtr problem = make(std::move(owned));
  const std::vector<double> base_x = base_point();
  std::vector<double> x(DIMENSION), y(1);

  if (!problem->FromBaseCoordinates(base_x, x))
    return 0;
  problem->EvaluateFunction(x, y);
  for (size_t i = 0; i < DIMENSION; ++i) {
    if (fabs(recorder->last_x[i] - base_x[i]) > precision * DoubleMax(1.0, fabs(base_x[i])))
      return 0;
  }
  return 1;
}

/** Checks that FromBaseCoordinates() reports that the problem built by make cannot be inverted. */
template <typename Make> static int is_not_invertible(Make make) {
  ProblemPtr problem = make(std::make_unique<Recorder>());
  const std::vector<double> base_x = base_point();
  std::vector<double> x(DIMENSION);
  return !problem->FromBaseCoordinates(base_x, x);
}

/**
 * Checks that the inverse of the transformation maps the image of some point back to the point, up to the given
 * relative precision (0 means exactly).
 */
template <typename Transformation> static int inverts(const Transformation &transformation, const double precision) {
  const std::vector<double> x = base_point();
  std::vector<double> y(DIMENSION), inverse_y(DIMENSION);
  transformation.Apply(x, y);
  transformation.Inverse(y, inverse_y);
  for (size_t i = 0; i < DIMENSION; ++i) {
    if (fabs(inverse_y[i] - x[i]) > precision * DoubleMax(1.0, fabs(x[i])))
      return 0;
  }
  return 1;
}

MU_TEST(test_exact_inverses) {
  const std::vector<double> offset = {0.5, -1.0, 2.0, -0.25, 3.0, 0.0};
  const size_t permutation[] = {3, 0, 5, 1, 4, 2};

  mu_check(maps_back([](ProblemPtr p) { return std::make_unique<TransformVarsEnforceBox>(std::move(p)); }, 0));
  mu_check(maps_back([](ProblemPtr p) { return std::make_unique<TransformObjScale>(std::move(p), 10.0); }, 0));

  mu_check(inverts(VarsShift(offset), 1e-15));
  mu_check(inverts(VarsScale(2.5, DIMENSION), 1e-15));
  mu_check(inverts(VarsConditioning(100.0, DIMENSION), 1e-15));
  mu_check(inverts(VarsPermutation(permutation, false), 0));
  mu_check(inverts(VarsPermutation(permutation, true), 0));
  mu_check(inverts(VarsXHat(4711, DIMENSION), 0));
  mu_check(inverts(VarsBrs(DIMENSION), 1e-15));
  mu_check(inverts(VarsZHat(offset), 1e-15));
}

MU_TEST(test_rotation_inverses) {
  size_t number_of_blocks;
  const std::vector<size_t> block_sizes = GetBlockSizes(&number_of_blocks, DIMENSION, "bbob-largescale");
  mu_check(inverts(VarsBlockrotation(ComputeBlockrotation(12345, block_sizes)), 1e-13));
  mu_check(inverts(VarsPermblockdiag(DIMENSION, 1, 2, 3), 1e-13));

  const Matrix R = Bbob2009ComputeRotation(12345, DIMENSION);
  const std::vector<double> b = {0.5, -1.0, 2.0, -0.25, 3.0, 0.0};
  mu_check(inverts(VarsAffine(R, b), 1e-13));

  /* A matrix that is not orthogonal is inverted by solving the linear system */
  Matrix M = R;
  M(0, 0) *= 2.0;
  mu_check(inverts(VarsAffine(M, b), 1e-12));

  /* A singular matrix is not invertible */
  for (size_t j = 0; j < DIMENSION; ++j)
    M(1, j) = M(0, j);
  std::vector<double> x(DIMENSION);
  mu_check(!VarsAffine(M, b).Inverse(base_point(), x));
}

MU_TEST(test_numerical_inverses) {
  mu_check(inverts(VarsOscillate(), 1e-9));
  mu_check(inverts(VarsAsymmetric(0.2, DIMENSION), 1e-9));
}

MU_TEST(test_chain_inverse) {
  mu_check(maps_back(
      [](ProblemPtr p) {
        p = std::make_unique<TransformObjScale>(std::move(p), 2.0);
        return std::make_unique<TransformVarsEnforceBox>(std::move(p));
      },
      0));
}

MU_TEST(test_rounding_is_not_invertible) {
  const std::vector<double> lower(DIMENSION, -5.0), upper(DIMENSION, 5.0);
  mu_check(is_not_invertible(
      [&](ProblemPtr p) { return std::make_unique<TransformVarsDiscretize>(std::move(p), lower, upper, DIMENSION); }));
}

} // namespace coco

/**
 * Run all tests in this file.
 */
int main(void) {
  MU_RUN_TEST(coco::test_exact_inverses);
  MU_RUN_TEST(coco::test_rotation_inverses);
  MU_RUN_TEST(coco::test_numerical_inverses);
  MU_RUN_TEST(coco::test_chain_inverse);
  MU_RUN_TEST(coco::test_rounding_is_not_invertible);

  MU_REPORT();

  return MU_EXIT_CODE;
}
