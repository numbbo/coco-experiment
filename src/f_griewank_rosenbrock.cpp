/**
 * @file f_griewank_rosenbrock.cpp
 * @brief Implementation of the Griewank-Rosenbrock function and problem.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

namespace coco {

/**
 * @brief The Griewank-Rosenbrock function.
 */
class GriewankRosenbrock final : public Problem {
public:
  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double Raw(std::span<const double> x) const;

  GriewankRosenbrock(const size_t number_of_variables, double facftrue) : Problem(number_of_variables, 1, 0) {
    InitFromScalars("Griewank Rosenbrock function", -5.0, 5.0, 1);
    SetId("%s_d%02zu", "griewank_rosenbrock", number_of_variables);

    this->facftrue_ = facftrue;

    /* Compute best solution */
    EvaluateFunction(best_parameter_, best_value_);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives_ == 1);
    y[0] = Raw(x);
    assert(y[0] + 1e-13 >= best_value_[0]);
  }

private:
  double facftrue_{};
};

/**
 * @brief Implements the Griewank-Rosenbrock function without connections to any
 * COCO structures.
 */
double GriewankRosenbrock::Raw(std::span<const double> x) const {

  size_t i = 0;
  double tmp = 0;

  /* Computation core */
  double result = 0.0;
  for (i = 0; i < number_of_variables_ - 1; ++i) {
    const double c1 = x[i] * x[i] - x[i + 1];
    const double c2 = 1.0 - x[i];
    tmp = 100.0 * c1 * c1 + c2 * c2;
    result += tmp / 4000. - cos(tmp);
  }
  result = facftrue_ + facftrue_ * result / static_cast<double>(number_of_variables_ - 1);

  return result;
}

/**
 * @brief The BBOB Griewank-Rosenbrock function f19, also in the bbob-boxed and bbob-noisy suites.
 *
 * f(x) = f_griewank_rosenbrock(s R x + 0.5) + fopt with s = max(1, sqrt(n) / 8).
 */
class BbobGriewankRosenbrock final : public BbobProblem<BbobGriewankRosenbrock> {
public:
  BbobGriewankRosenbrock(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                         const double facftrue, const char *problem_id_template, const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "4-multi-modal", problem_id_template, problem_name_template),
        rotation_(ScaledRotationMatrix(rseed, dimension), std::vector<double>(dimension)),
        shift_(std::vector<double>(dimension, -0.5)), griewank_rosenbrock_(dimension, facftrue),
        fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension), z2_(dimension) {
    /* The optimum of the Griewank-Rosenbrock function is (1, ..., 1), map it back through the transformations */
    const std::vector<double> ones(dimension, 1.0);
    shift_.Inverse(ones, z1_);
    rotation_.Inverse(z1_, best_parameter_);
    /* Evaluating at the best parameter could add rounding errors */
    best_value_[0] = griewank_rosenbrock_.Raw(ones) + fopt_;
  }

  double operator()(std::span<const double> x) {
    rotation_.Apply(x, z1_);
    shift_.Apply(z1_, z2_);
    return griewank_rosenbrock_.Raw(z2_) + fopt_;
  }

private:
  /** @brief Returns s R with the rotation R. */
  static Matrix ScaledRotationMatrix(const long rseed, const size_t dimension) {
    Matrix M = Bbob2009ComputeRotation(rseed, dimension);
    const double scales = DoubleMax(1., sqrt(static_cast<double>(dimension)) / 8.);
    for (size_t i = 0; i < dimension; ++i) {
      for (size_t j = 0; j < dimension; ++j)
        M(i, j) *= scales;
    }
    return M;
  }

  const VarsAffine rotation_;
  const VarsShift shift_; /**< @brief Adds 0.5 */
  const GriewankRosenbrock griewank_rosenbrock_;
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

/**
 * @brief The Griewank-Rosenbrock function f19 of the bbob-largescale suite, with a permuted block rotation.
 *
 * f(x) = f_griewank_rosenbrock(s R x + 0.5) + fopt with the permuted block rotation R and
 * s = max(1, sqrt(block size) / 8).
 */
class BbobGriewankRosenbrockPermblockdiag final : public BbobProblem<BbobGriewankRosenbrockPermblockdiag> {
public:
  BbobGriewankRosenbrockPermblockdiag(const size_t function, const size_t dimension, const size_t instance,
                                      const long rseed, const double facftrue, const char *problem_id_template,
                                      const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "4-multi-modal", problem_id_template, problem_name_template),
        rotation_(dimension, rseed + 2000000, rseed, rseed + 3000000),
        scale_(DoubleMax(1.0, sqrt(static_cast<double>(RotationMatrixBlockSize(dimension))) / 8.0), dimension),
        shift_(std::vector<double>(dimension, -0.5)), griewank_rosenbrock_(dimension, facftrue),
        fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension), z2_(dimension) {
    /* The optimum of the Griewank-Rosenbrock function is (1, ..., 1), map it back through the transformations */
    const std::vector<double> ones(dimension, 1.0);
    shift_.Inverse(ones, z1_);
    scale_.Inverse(z1_, z1_);
    rotation_.Inverse(z1_, best_parameter_);
    /* Evaluating at the best parameter could add rounding errors */
    best_value_[0] = griewank_rosenbrock_.Raw(ones) + fopt_;
  }

  double operator()(std::span<const double> x) {
    rotation_.Apply(x, z1_);
    scale_.Apply(z1_, z1_);
    shift_.Apply(z1_, z2_);
    return griewank_rosenbrock_.Raw(z2_) + fopt_;
  }

private:
  const VarsPermblockdiag rotation_;
  const VarsScale scale_;
  const VarsShift shift_; /**< @brief Adds 0.5 */
  const GriewankRosenbrock griewank_rosenbrock_;
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

} // namespace coco
