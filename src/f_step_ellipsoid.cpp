/**
 * @file f_step_ellipsoid.cpp
 * @brief Implementation of the step ellipsoid function and problem.
 *
 * The BBOB step ellipsoid function intertwines the variable and objective
 * transformations in such a way that it is hard to devise a composition of
 * generic transformations to implement it. In the end one would have to
 * implement several custom transformations which would be used solely by this
 * problem. Therefore we opt to implement it as a monolithic function instead.
 *
 * TODO: It would be nice to have a generic step ellipsoid function to
 * complement this one.
 */
#include <assert.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "coco_utilities.cpp"
#include "suite_bbob_legacy_code.cpp"

#include "bbob_problem.cpp"
#include "transform_functions.cpp"

namespace coco {

/**
 * @brief The BBOB step ellipsoid function f7, also in the bbob-boxed and bbob-noisy suites.
 *
 * f(x) = 0.1 max(|ẑ_1| 1e-4, f_ellipsoid^100(Q round(ẑ))) + penalty_scale f_pen(x) + fopt with
 * ẑ = Λ^10 R (x - xopt) and the rounding to steps of 1 or 0.1 (for |ẑ_i| <= 0.5).
 */
class BbobStepEllipsoid final : public BbobProblem<BbobStepEllipsoid> {
public:
  BbobStepEllipsoid(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                    const BbobVariant variant, const double penalty_scale, const char *problem_id_template,
                    const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "2-moderate", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)),
        conditioned_rotation_(ConditionedRotationMatrix(rseed, dimension), std::vector<double>(dimension)),
        round_step_(10.0),
        rotation_(Bbob2009ComputeRotation(rseed + 1000000, dimension), std::vector<double>(dimension)),
        scales_(dimension), fopt_(Bbob2009ComputeFopt(function, instance)), penalty_scale_(penalty_scale),
        z1_(dimension), z2_(dimension) {
    static const double condition = 100;
    assert(dimension > 1);
    for (size_t i = 0; i < dimension; ++i)
      scales_[i] = pow(condition, static_cast<double>(i) / (static_cast<double>(dimension) - 1.0));

    /* The optimum is at xopt, where ẑ = 0 */
    shift_.Inverse(best_parameter_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty =
        penalty_scale_ * BoundaryPenalty(x, smallest_values_of_interest_, largest_values_of_interest_);
    shift_.Apply(x, z1_);
    conditioned_rotation_.Apply(z1_, z2_);
    const double zhat_1 = z2_[0];
    round_step_.Apply(z2_, z2_);
    rotation_.Apply(z2_, z1_);

    double y = 0.0;
    for (size_t i = 0; i < number_of_variables_; ++i)
      y += scales_[i] * z1_[i] * z1_[i];
    return 0.1 * DoubleMax(fabs(zhat_1) * 1.0e-4, y) + penalty + fopt_;
  }

private:
  /** @brief Returns Λ^10 R with the rotation R. */
  static Matrix ConditionedRotationMatrix(const long rseed, const size_t dimension) {
    static const double condition = 100;
    Matrix M = Bbob2009ComputeRotation(rseed, dimension);
    for (size_t i = 0; i < dimension; ++i) {
      const double c1 = sqrt(pow(condition / 10., static_cast<double>(i) / static_cast<double>(dimension - 1)));
      for (size_t j = 0; j < dimension; ++j)
        M(i, j) = c1 * M(i, j);
    }
    return M;
  }

  const VarsShift shift_;
  const VarsAffine conditioned_rotation_;
  const VarsRoundStep round_step_;
  const VarsAffine rotation_;
  std::vector<double> scales_; /**< @brief condition^(i / (n - 1)) */
  const double fopt_;
  const double penalty_scale_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

/**
 * @brief The step ellipsoid function f7 of the bbob-largescale suite, with permuted block rotations.
 *
 * f(x) = 0.1 max(|ẑ_1| n 1e-4, f_ellipsoid^100(Q round(ẑ))) / n + f_pen(x) + fopt with ẑ = Λ^10 R (x - xopt), the
 * rounding to steps of 1 or 0.1 (for |ẑ_i| <= 0.5), and the permuted block rotations R and Q.
 */
class BbobStepEllipsoidPermblockdiag final : public BbobProblem<BbobStepEllipsoidPermblockdiag> {
public:
  BbobStepEllipsoidPermblockdiag(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                                 const BbobVariant variant, const char *problem_id_template,
                                 const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "2-moderate", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)), rotation1_(dimension, rseed + 2000000, rseed, rseed + 3000000),
        conditioning_(10.0, dimension), round_step_(10.0),
        rotation2_(dimension, rseed + 4000000, rseed + 1000000, rseed + 5000000), scales_(dimension),
        fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension), z2_(dimension) {
    static const double condition = 100;
    for (size_t i = 0; i < dimension; ++i)
      scales_[i] = pow(condition, static_cast<double>(i) / (static_cast<double>(dimension) - 1.0));

    /* The optimum of the ellipsoid is the origin, map it back through the transformations; the rounding keeps it */
    rotation2_.Inverse(best_parameter_, z1_);
    round_step_.Apply(z1_, z1_);
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
    /* multiplied by the dimension to counter-balance the normalization */
    const double zhat_1 = fabs(z2_[0]) * static_cast<double>(number_of_variables_);
    round_step_.Apply(z2_, z2_);
    rotation2_.Apply(z2_, z1_);

    double y = 0.0;
    for (size_t i = 0; i < number_of_variables_; ++i)
      y += scales_[i] * z1_[i] * z1_[i];
    y = 0.1 * DoubleMax(zhat_1 * 1.0e-4, y);
    return y * NormByDimFactor(number_of_variables_) + PENALTY_FACTOR * penalty + fopt_;
  }

private:
  static constexpr double PENALTY_FACTOR = 1.0;

  const VarsShift shift_;
  const VarsPermblockdiag rotation1_;
  const VarsConditioning conditioning_;
  const VarsRoundStep round_step_;
  const VarsPermblockdiag rotation2_;
  std::vector<double> scales_; /**< @brief condition^(i / (n - 1)) */
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

} // namespace coco
