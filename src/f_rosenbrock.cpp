/**
 * @file f_rosenbrock.cpp
 * @brief Implementation of the Rosenbrock function and problem.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

namespace coco {

/**
 * @brief Implements the Rosenbrock function without connections to any COCO
 * structures.
 */
static double RosenbrockRaw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  size_t i = 0;
  double s1 = 0.0, s2 = 0.0, tmp;

  assert(number_of_variables > 1);

  for (i = 0; i < number_of_variables - 1; ++i) {
    tmp = (x[i] * x[i] - x[i + 1]);
    s1 += tmp * tmp;
    tmp = (x[i] - 1.0);
    s2 += tmp * tmp;
  }
  double result = 100.0 * s1 + s2;

  return result;
}

/**
 * @brief The Rosenbrock function.
 */
class Rosenbrock final : public Problem {
public:
  explicit Rosenbrock(const size_t number_of_variables) : Problem(number_of_variables, 1, 0) {
    InitFromScalars("Rosenbrock function", -5.0, 5.0, 1.0);
    SetId("%s_d%02zu", "rosenbrock", number_of_variables);

    /* Compute best solution */
    EvaluateFunction(best_parameter_, best_value_);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives_ == 1);
    y[0] = RosenbrockRaw(x);
    assert(y[0] + 1e-13 >= best_value_[0]);
  }
};

/**
 * @brief The location of the optimum of the BBOB Rosenbrock functions before the shift to (1, ..., 1): 0.75 xopt.
 */
static std::vector<double> RosenbrockXopt(const long rseed, const size_t dimension, const BbobVariant variant) {
  std::vector<double> xopt = BbobXopt(rseed, dimension, variant);
  for (size_t i = 0; i < dimension; ++i)
    xopt[i] *= 0.75;
  return xopt;
}

/**
 * @brief The scaling factor of the BBOB Rosenbrock functions: max(1, sqrt(n) / 8), where n is the block size in
 * the large-scale suite.
 */
static double RosenbrockFactor(const size_t dimension, const BbobVariant variant) {
  if (variant == BbobVariant::LARGESCALE) {
    const size_t block_size = RotationMatrixBlockSize(dimension);
    return DoubleMax(1.0, sqrt(static_cast<double>(block_size)) / 8.0);
  }
  return DoubleMax(1.0, sqrt(static_cast<double>(dimension)) / 8.0);
}

/**
 * @brief The BBOB Rosenbrock function f8, also in the bbob-boxed, bbob-largescale and bbob-noisy suites.
 *
 * f(x) = f_rosenbrock(s (x - 0.75 xopt) + 1) + fopt with s = max(1, sqrt(n) / 8), where n is the block size in the
 * large-scale suite; normalized by the dimension before fopt is added in the large-scale suite.
 */
class BbobRosenbrock final : public BbobProblem<BbobRosenbrock> {
public:
  BbobRosenbrock(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                 const BbobVariant variant, const char *problem_id_template, const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "2-moderate", problem_id_template, problem_name_template),
        shift_(RosenbrockXopt(rseed, dimension, variant)), scale_(RosenbrockFactor(dimension, variant), dimension),
        shift_to_one_(std::vector<double>(dimension, -1.0)), fopt_(Bbob2009ComputeFopt(function, instance)),
        normalize_(variant == BbobVariant::LARGESCALE), z_(dimension) {
    /* The optimum of the Rosenbrock function is (1, ..., 1), map it back through the transformations */
    std::fill(z_.begin(), z_.end(), 1.0);
    shift_to_one_.Inverse(z_, z_);
    scale_.Inverse(z_, z_);
    shift_.Inverse(z_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    shift_.Apply(x, z_);
    scale_.Apply(z_, z_);
    shift_to_one_.Apply(z_, z_);
    double y = RosenbrockRaw(z_);
    if (normalize_)
      y *= NormByDimFactor(number_of_variables_);
    return y + fopt_;
  }

private:
  const VarsShift shift_;
  const VarsScale scale_;
  const VarsShift shift_to_one_; /**< @brief Adds 1 */
  const double fopt_;
  const bool normalize_;
  std::vector<double> z_; /**< @brief The transformed decision values */
};

/**
 * @brief The BBOB rotated Rosenbrock function f9, also in the bbob-boxed suite.
 *
 * f(x) = f_rosenbrock(s R x + 0.5) + fopt with s = max(1, sqrt(n) / 8).
 */
class BbobRosenbrockRotated final : public BbobProblem<BbobRosenbrockRotated> {
public:
  BbobRosenbrockRotated(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                        const char *problem_id_template, const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "2-moderate", problem_id_template, problem_name_template),
        affine_(ScaledRotationMatrix(rseed, dimension), std::vector<double>(dimension, 0.5)),
        fopt_(Bbob2009ComputeFopt(function, instance)), z_(dimension) {
    /* The optimum of the Rosenbrock function is (1, ..., 1), map it back through the transformation */
    const std::vector<double> ones(dimension, 1.0);
    affine_.Inverse(ones, best_parameter_);
    /* Evaluating at the best parameter could add rounding errors */
    best_value_[0] = RosenbrockRaw(ones) + fopt_;
  }

  double operator()(std::span<const double> x) {
    affine_.Apply(x, z_);
    return RosenbrockRaw(z_) + fopt_;
  }

private:
  /** @brief Returns s R with the rotation R. */
  static Matrix ScaledRotationMatrix(const long rseed, const size_t dimension) {
    const double factor = DoubleMax(1.0, sqrt(static_cast<double>(dimension)) / 8.0);
    Matrix M = Bbob2009ComputeRotation(rseed, dimension);
    for (size_t row = 0; row < dimension; ++row) {
      for (size_t column = 0; column < dimension; ++column)
        M(row, column) = factor * M(row, column);
    }
    return M;
  }

  const VarsAffine affine_;
  const double fopt_;
  std::vector<double> z_; /**< @brief The transformed decision values */
};

/**
 * @brief The rotated Rosenbrock function f9 of the bbob-largescale suite, with a permuted block rotation.
 *
 * f(x) = f_rosenbrock(s R (x - 0.75 xopt) + 1) / n + fopt with the permuted block rotation R and
 * s = max(1, sqrt(block size) / 8).
 */
class BbobRosenbrockPermblockdiag final : public BbobProblem<BbobRosenbrockPermblockdiag> {
public:
  BbobRosenbrockPermblockdiag(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                              const BbobVariant variant, const char *problem_id_template,
                              const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "2-moderate", problem_id_template, problem_name_template),
        shift_(RosenbrockXopt(rseed, dimension, variant)),
        rotation_(dimension, rseed + 2000000, rseed, rseed + 3000000),
        scale_(RosenbrockFactor(dimension, BbobVariant::LARGESCALE), dimension),
        shift_to_one_(std::vector<double>(dimension, -1.0)), fopt_(Bbob2009ComputeFopt(function, instance)),
        z1_(dimension), z2_(dimension) {
    /* The optimum of the Rosenbrock function is (1, ..., 1), map it back through the transformations */
    std::fill(z1_.begin(), z1_.end(), 1.0);
    shift_to_one_.Inverse(z1_, z1_);
    scale_.Inverse(z1_, z1_);
    rotation_.Inverse(z1_, z2_);
    shift_.Inverse(z2_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    shift_.Apply(x, z1_);
    rotation_.Apply(z1_, z2_);
    scale_.Apply(z2_, z2_);
    shift_to_one_.Apply(z2_, z2_);
    return RosenbrockRaw(z2_) * NormByDimFactor(number_of_variables_) + fopt_;
  }

private:
  const VarsShift shift_;
  const VarsPermblockdiag rotation_;
  const VarsScale scale_;
  const VarsShift shift_to_one_; /**< @brief Adds 1 */
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

} // namespace coco
