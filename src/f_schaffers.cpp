/**
 * @file f_schaffers.cpp
 * @brief Implementation of the Schaffer's F7 function and problem,
 * transformations not implemented for the moment.
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
 * @brief Implements the Schaffer's F7 function without connections to any COCO
 * structures.
 */
static double SchaffersRaw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  size_t i = 0;

  assert(number_of_variables > 1);

  /* Computation core */
  double result = 0.0;
  for (i = 0; i < number_of_variables - 1; ++i) {
    const double tmp = x[i] * x[i] + x[i + 1] * x[i + 1];
    if (IsInf(tmp) && IsNan(sin(50.0 * pow(tmp, 0.1)))) /* sin(inf) -> nan */
      /* the second condition is necessary to pass the integration tests under Windows and Linux */
      return tmp;
    result += pow(tmp, 0.25) * (1.0 + pow(sin(50.0 * pow(tmp, 0.1)), 2.0));
  }
  result = pow(result / (static_cast<double>(number_of_variables) - 1.0), 2.0);

  return result;
}

/**
 * @brief The BBOB Schaffer's F7 functions f17 and f18, also in the bbob-boxed and bbob-noisy suites.
 *
 * f(x) = f_schaffers(Λ^conditioning Q T_asy^0.5(R (x - xopt))) + fopt + penalty_scale f_pen(x).
 */
class BbobSchaffers final : public BbobProblem<BbobSchaffers> {
public:
  BbobSchaffers(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                const BbobVariant variant, const double conditioning, const double penalty_scale,
                const char *problem_id_template, const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "4-multi-modal", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)),
        rotation_(Bbob2009ComputeRotation(rseed + 1000000, dimension), std::vector<double>(dimension)),
        asymmetric_(0.5, dimension), conditioned_rotation_(ConditionedRotationMatrix(rseed, dimension, conditioning),
                                                           std::vector<double>(dimension)),
        fopt_(Bbob2009ComputeFopt(function, instance)), penalty_scale_(penalty_scale), z1_(dimension), z2_(dimension) {
    /* The optimum of Schaffer's function is the origin, map it back through the transformations */
    conditioned_rotation_.Inverse(best_parameter_, z1_);
    asymmetric_.Inverse(z1_, z1_);
    rotation_.Inverse(z1_, z2_);
    shift_.Inverse(z2_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = BoundaryPenalty(x, smallest_values_of_interest_, largest_values_of_interest_);
    shift_.Apply(x, z1_);
    rotation_.Apply(z1_, z2_);
    asymmetric_.Apply(z2_, z2_);
    conditioned_rotation_.Apply(z2_, z1_);
    const double y = SchaffersRaw(z1_) + fopt_;
    return y + penalty_scale_ * penalty;
  }

private:
  /** @brief Returns Λ^conditioning Q with the rotation Q. */
  static Matrix ConditionedRotationMatrix(const long rseed, const size_t dimension, const double conditioning) {
    Matrix M = Bbob2009ComputeRotation(rseed, dimension);
    for (size_t i = 0; i < dimension; ++i) {
      for (size_t j = 0; j < dimension; ++j) {
        double exponent = 1.0 * static_cast<int>(i) / (static_cast<double>(dimension) - 1.0);
        M(i, j) = M(i, j) * pow(sqrt(conditioning), exponent);
      }
    }
    return M;
  }

  const VarsShift shift_;
  const VarsAffine rotation_;
  const VarsAsymmetric asymmetric_;
  const VarsAffine conditioned_rotation_;
  const double fopt_;
  const double penalty_scale_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

/**
 * @brief The Schaffer's F7 functions f17 and f18 of the bbob-largescale suite, with permuted block rotations.
 *
 * f(x) = f_schaffers(Λ^conditioning Q T_asy^0.5(R (x - xopt))) + 10 f_pen(x) + fopt with the permuted block rotations
 * R and Q.
 */
class BbobSchaffersPermblockdiag final : public BbobProblem<BbobSchaffersPermblockdiag> {
public:
  BbobSchaffersPermblockdiag(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                             const BbobVariant variant, const double condition, const char *problem_id_template,
                             const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "4-multi-modal", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)),
        rotation1_(dimension, rseed + 5000000, rseed + 1000000, rseed + 6000000), asymmetric_(0.5, dimension),
        rotation2_(dimension, rseed + 3000000, rseed, rseed + 4000000), conditioning_(condition, dimension),
        fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension), z2_(dimension) {
    /* The optimum of Schaffer's function is the origin, map it back through the transformations */
    conditioning_.Inverse(best_parameter_, z1_);
    rotation2_.Inverse(z1_, z2_);
    asymmetric_.Inverse(z2_, z2_);
    rotation1_.Inverse(z2_, z1_);
    shift_.Inverse(z1_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = BoundaryPenalty(x, smallest_values_of_interest_, largest_values_of_interest_);
    shift_.Apply(x, z1_);
    rotation1_.Apply(z1_, z2_);
    asymmetric_.Apply(z2_, z2_);
    rotation2_.Apply(z2_, z1_);
    conditioning_.Apply(z1_, z1_);
    const double y = SchaffersRaw(z1_) + PENALTY_FACTOR * penalty;
    return y + fopt_;
  }

private:
  static constexpr double PENALTY_FACTOR = 10.0;

  const VarsShift shift_;
  const VarsPermblockdiag rotation1_;
  const VarsAsymmetric asymmetric_;
  const VarsPermblockdiag rotation2_;
  const VarsConditioning conditioning_;
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

} // namespace coco
