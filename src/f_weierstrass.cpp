/**
 * @file f_weierstrass.cpp
 * @brief Implementation of the Weierstrass function and problem.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

namespace coco {

/** @brief Number of summands in the Weierstrass problem. */
static constexpr size_t F_WEIERSTRASS_SUMMANDS = 12;

/**
 * @brief The Weierstrass function.
 */
class Weierstrass final : public Problem {
public:
  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double Raw(std::span<const double> x) const;

  explicit Weierstrass(const size_t number_of_variables) : Problem(number_of_variables, 1, 0) {

    InitFromScalars("Weierstrass function", -5.0, 5.0, 0.);
    SetId("%s_d%02zu", "weierstrass", number_of_variables);

    f0_ = 0.0;
    for (size_t i = 0; i < F_WEIERSTRASS_SUMMANDS; ++i) {
      ak_[i] = pow(0.5, static_cast<double>(i));
      bk_[i] = pow(3., static_cast<double>(i));
      f0_ += ak_[i] * cos(2 * coco_pi * bk_[i] * 0.5);
    }

    /* Compute best solution */
    std::vector<double> non_unique_best_value(number_of_variables, 0.0);
    EvaluateFunction(non_unique_best_value, best_value_);
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
  double f0_{};
  double ak_[F_WEIERSTRASS_SUMMANDS]{};
  double bk_[F_WEIERSTRASS_SUMMANDS]{};
};

/**
 * @brief Implements the Weierstrass function without connections to any COCO
 * structures.
 */
double Weierstrass::Raw(std::span<const double> x) const {

  double result = 0.0;
  for (size_t i = 0; i < number_of_variables_; ++i) {
    for (size_t j = 0; j < F_WEIERSTRASS_SUMMANDS; ++j) {
      result += cos(2 * coco_pi * (x[i] + 0.5) * bk_[j]) * ak_[j];
    }
  }
  result = 10.0 * pow(result / static_cast<double>(number_of_variables_) - f0_, 3.0);

  return result;
}

/**
 * @brief The BBOB Weierstrass function f16, also in the bbob-boxed suite.
 *
 * f(x) = f_weierstrass(R Λ^(1/100) Q T_osz(R (x - xopt))) + fopt + 10 / n f_pen(x).
 */
class BbobWeierstrass final : public BbobProblem<BbobWeierstrass> {
public:
  BbobWeierstrass(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                  const BbobVariant variant, const char *problem_id_template, const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "4-multi-modal", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)),
        rotation_(Bbob2009ComputeRotation(rseed + 1000000, dimension), std::vector<double>(dimension)),
        conditioned_rotation_(ConditionedRotationMatrix(rseed, dimension), std::vector<double>(dimension)),
        weierstrass_(dimension), fopt_(Bbob2009ComputeFopt(function, instance)),
        penalty_factor_(10.0 / static_cast<double>(dimension)), z1_(dimension), z2_(dimension) {
    /* The optimum of the Weierstrass function is the origin, map it back through the transformations */
    conditioned_rotation_.Inverse(best_parameter_, z1_);
    oscillate_.Inverse(z1_, z1_);
    rotation_.Inverse(z1_, z2_);
    shift_.Inverse(z2_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = BoundaryPenalty(x, smallest_values_of_interest_, largest_values_of_interest_);
    shift_.Apply(x, z1_);
    rotation_.Apply(z1_, z2_);
    oscillate_.Apply(z2_, z2_);
    conditioned_rotation_.Apply(z2_, z1_);
    const double y = weierstrass_.Raw(z1_) + fopt_;
    return y + penalty_factor_ * penalty;
  }

private:
  /** @brief Returns R Λ^(1/100) Q with the rotations R and Q. */
  static Matrix ConditionedRotationMatrix(const long rseed, const size_t dimension) {
    const double condition = 100.0;
    const Matrix rot1 = Bbob2009ComputeRotation(rseed + 1000000, dimension);
    const Matrix rot2 = Bbob2009ComputeRotation(rseed, dimension);
    std::vector<double> scales(dimension);
    for (size_t k = 0; k < dimension; ++k) {
      const double base = 1.0 / sqrt(condition);
      const double exponent = 1.0 * static_cast<int>(k) / (static_cast<double>(dimension) - 1.0);
      scales[k] = pow(base, exponent);
    }
    return MatrixProductWithDiagonal(rot1, scales, rot2);
  }

  const VarsShift shift_;
  const VarsAffine rotation_;
  const VarsOscillate oscillate_{};
  const VarsAffine conditioned_rotation_;
  const Weierstrass weierstrass_;
  const double fopt_;
  const double penalty_factor_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

/**
 * @brief The Weierstrass function f16 of the bbob-largescale suite, with permuted block rotations.
 *
 * f(x) = f_weierstrass(R Λ^(1/100) Q T_osz(R (x - xopt))) + 10 / n f_pen(x) + fopt with the permuted block rotations
 * R and Q.
 */
class BbobWeierstrassPermblockdiag final : public BbobProblem<BbobWeierstrassPermblockdiag> {
public:
  BbobWeierstrassPermblockdiag(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                               const BbobVariant variant, const char *problem_id_template,
                               const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "4-multi-modal", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)),
        rotation1_(dimension, rseed + 5000000, rseed + 1000000, rseed + 6000000),
        rotation2_(dimension, rseed + 3000000, rseed, rseed + 4000000), conditioning_(1.0 / 100.0, dimension),
        weierstrass_(dimension), fopt_(Bbob2009ComputeFopt(function, instance)),
        penalty_factor_(10.0 / static_cast<double>(dimension)), z1_(dimension), z2_(dimension) {
    /* The optimum of the Weierstrass function is the origin, map it back through the transformations */
    rotation1_.Inverse(best_parameter_, z1_);
    conditioning_.Inverse(z1_, z1_);
    rotation2_.Inverse(z1_, z2_);
    oscillate_.Inverse(z2_, z2_);
    rotation1_.Inverse(z2_, z1_);
    shift_.Inverse(z1_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = BoundaryPenalty(x, smallest_values_of_interest_, largest_values_of_interest_);
    shift_.Apply(x, z1_);
    rotation1_.Apply(z1_, z2_);
    oscillate_.Apply(z2_, z2_);
    rotation2_.Apply(z2_, z1_);
    conditioning_.Apply(z1_, z1_);
    rotation1_.Apply(z1_, z2_);
    const double y = weierstrass_.Raw(z2_) + penalty_factor_ * penalty;
    return y + fopt_;
  }

private:
  const VarsShift shift_;
  const VarsPermblockdiag rotation1_;
  const VarsOscillate oscillate_{};
  const VarsPermblockdiag rotation2_;
  const VarsConditioning conditioning_;
  const Weierstrass weierstrass_;
  const double fopt_;
  const double penalty_factor_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

} // namespace coco
