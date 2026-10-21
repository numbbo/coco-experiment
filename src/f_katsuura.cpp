/**
 * @file f_katsuura.cpp
 * @brief Implementation of the Katsuura function and problem.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "coco_utilities.cpp"
#include "suite_bbob_legacy_code.cpp"

namespace coco {

/**
 * @brief Implements the Katsuura function without connections to any COCO
 * structures.
 */
static double KatsuuraRaw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  const double exponent = 10. / pow(static_cast<double>(number_of_variables), 1.2);

  /* Computation core */
  double result = 1.0;
  for (size_t i = 0; i < number_of_variables; ++i) {
    double tmp = 0;
    for (int j = 1; j < 33; ++j) {
      const double tmp2 = ldexp(1.0, j); /* 2^j */
      tmp += fabs(tmp2 * x[i] - DoubleRound(tmp2 * x[i])) / tmp2;
    }
    tmp = 1.0 + (static_cast<double>(i) + 1) * tmp;
    result *= pow(tmp, exponent);
  }
  result =
      10. / (static_cast<double>(number_of_variables)) / (static_cast<double>(number_of_variables)) * (-1. + result);

  return result;
}

/**
 * @brief The BBOB Katsuura function f23, also in the bbob-boxed suite.
 *
 * f(x) = f_katsuura(Q Λ^100 R (x - xopt)) + fopt + f_pen(x).
 */
class BbobKatsuura final : public BbobProblem<BbobKatsuura> {
public:
  BbobKatsuura(const size_t function, const size_t dimension, const size_t instance, const long rseed,
               const BbobVariant variant, const char *problem_id_template, const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "5-weakly-structured", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)),
        conditioned_rotation_(ConditionedRotationMatrix(rseed, dimension), std::vector<double>(dimension)),
        fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension), z2_(dimension) {
    /* The optimum of the Katsuura function is the origin, map it back through the transformations */
    conditioned_rotation_.Inverse(best_parameter_, z1_);
    shift_.Inverse(z1_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = BoundaryPenalty(x, smallest_values_of_interest_, largest_values_of_interest_);
    shift_.Apply(x, z1_);
    conditioned_rotation_.Apply(z1_, z2_);
    /* There is no shift 'fopt' in the definition */
    const double y = KatsuuraRaw(z2_) + fopt_;
    return y + PENALTY_FACTOR * penalty;
  }

private:
  /** @brief Returns Q Λ^100 R with the rotations Q and R. */
  static Matrix ConditionedRotationMatrix(const long rseed, const size_t dimension) {
    const Matrix rot1 = Bbob2009ComputeRotation(rseed + 1000000, dimension);
    const Matrix rot2 = Bbob2009ComputeRotation(rseed, dimension);
    std::vector<double> scales(dimension);
    for (size_t k = 0; k < dimension; ++k) {
      const double exponent = 1.0 * static_cast<int>(k) / (static_cast<double>(dimension) - 1.0);
      scales[k] = pow(sqrt(100), exponent);
    }
    return MatrixProductWithDiagonal(rot1, scales, rot2);
  }

  static constexpr double PENALTY_FACTOR = 1.0;

  const VarsShift shift_;
  const VarsAffine conditioned_rotation_;
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

/**
 * @brief The Katsuura function f23 of the bbob-largescale suite, with permuted block rotations.
 *
 * f(x) = f_katsuura(Q Λ^100 R (x - xopt)) + f_pen(x) + fopt with the permuted block rotations R and Q.
 */
class BbobKatsuuraPermblockdiag final : public BbobProblem<BbobKatsuuraPermblockdiag> {
public:
  BbobKatsuuraPermblockdiag(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                            const BbobVariant variant, const char *problem_id_template,
                            const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "5-weakly-structured", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)), rotation1_(dimension, rseed + 3000000, rseed, rseed + 4000000),
        conditioning_(100.0, dimension), rotation2_(dimension, rseed + 5000000, rseed + 1000000, rseed + 6000000),
        fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension), z2_(dimension) {
    /* The optimum of the Katsuura function is the origin, map it back through the transformations */
    rotation2_.Inverse(best_parameter_, z1_);
    conditioning_.Inverse(z1_, z1_);
    rotation1_.Inverse(z1_, z2_);
    shift_.Inverse(z2_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = BoundaryPenalty(x, smallest_values_of_interest_, largest_values_of_interest_);
    shift_.Apply(x, z1_);
    rotation1_.Apply(z1_, z2_);
    conditioning_.Apply(z2_, z2_);
    rotation2_.Apply(z2_, z1_);
    const double y = KatsuuraRaw(z1_) + PENALTY_FACTOR * penalty;
    /* TODO: documentation, there is no fopt in the definition of this function */
    return y + fopt_;
  }

private:
  static constexpr double PENALTY_FACTOR = 1.0;

  const VarsShift shift_;
  const VarsPermblockdiag rotation1_;
  const VarsConditioning conditioning_;
  const VarsPermblockdiag rotation2_;
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

} // namespace coco
