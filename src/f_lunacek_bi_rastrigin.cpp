/**
 * @file f_lunacek_bi_rastrigin.cpp
 * @brief Implementation of the Lunacek bi-Rastrigin function and problem.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "bbob_problem.cpp"
#include "f_sphere.cpp"

namespace coco {

/**
 * @brief The BBOB Lunacek bi-Rastrigin function f24, also in the bbob-boxed suite.
 *
 * f(x) = min(Σ (x^_i - μ0)^2, d n + s Σ (x^_i - μ1)^2) + 10 (n - Σ cos(2π z_i)) + 10^4 f_pen(x) + fopt with
 * x^ = 2 sign(xopt) x and z = Q Λ^100 R (x^ - μ0).
 */
class BbobLunacekBiRastrigin final : public BbobProblem<BbobLunacekBiRastrigin> {
public:
  BbobLunacekBiRastrigin(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                         const BbobVariant variant, const char *problem_id_template, const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "5-weakly-structured", problem_id_template, problem_name_template),
        xopt_(LunacekXopt(rseed, dimension)), shift_(std::vector<double>(dimension, MU0)),
        conditioned_rotation_(ConditionedRotationMatrix(rseed, dimension), std::vector<double>(dimension)),
        rotation_(Bbob2009ComputeRotation(rseed + 1000000, dimension), std::vector<double>(dimension)),
        s_(1. - 0.5 / (sqrt(static_cast<double>(dimension + 20)) - 4.1)), mu1_(-sqrt((MU0 * MU0 - D) / s_)),
        fopt_(Bbob2009ComputeFopt(function, instance)), x_hat_(dimension), z1_(dimension), z2_(dimension) {
    /* The variant only matters for the xopt that the legacy code computes and discards */
    (void)variant;
    best_parameter_ = xopt_;
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = BoundaryPenalty(x, smallest_values_of_interest_, largest_values_of_interest_);
    double sum1 = 0., sum2 = 0., sum3 = 0.;

    for (size_t i = 0; i < number_of_variables_; ++i) {
      x_hat_[i] = 2. * x[i];
      if (xopt_[i] < 0.) {
        x_hat_[i] *= -1.;
      }
    }
    shift_.Apply(x_hat_, z1_);
    conditioned_rotation_.Apply(z1_, z2_);
    rotation_.Apply(z2_, z1_);
    for (size_t i = 0; i < number_of_variables_; ++i) {
      sum1 += (x_hat_[i] - MU0) * (x_hat_[i] - MU0);
      sum2 += (x_hat_[i] - mu1_) * (x_hat_[i] - mu1_);
      sum3 += cos(2 * coco_pi * z1_[i]);
    }
    const double y = DoubleMin(sum1, D * static_cast<double>(number_of_variables_) + s_ * sum2) +
                     10. * (static_cast<double>(number_of_variables_) - sum3) + 1e4 * penalty;
    return y + fopt_;
  }

private:
  static constexpr double MU0 = 2.5;
  static constexpr double D = 1.;

  /** @brief Returns xopt = ±μ0 / 2 with random signs. */
  static std::vector<double> LunacekXopt(const long rseed, const size_t dimension) {
    std::vector<double> signs(dimension), xopt(dimension);
    Bbob2009Gauss(signs, rseed);
    for (size_t i = 0; i < dimension; ++i) {
      xopt[i] = 0.5 * MU0;
      if (signs[i] < 0.0) {
        xopt[i] *= -1.0;
      }
    }
    return xopt;
  }

  /** @brief Returns Λ^100 R with the rotation R. */
  static Matrix ConditionedRotationMatrix(const long rseed, const size_t dimension) {
    static const double condition = 100.;
    assert(dimension > 1);
    Matrix M = Bbob2009ComputeRotation(rseed, dimension);
    for (size_t i = 0; i < dimension; ++i) {
      const double c1 = pow(sqrt(condition), (static_cast<double>(i)) / static_cast<double>(dimension - 1));
      for (size_t j = 0; j < dimension; ++j)
        M(i, j) = c1 * M(i, j);
    }
    return M;
  }

  const std::vector<double> xopt_;
  const VarsShift shift_; /**< @brief Subtracts μ0 */
  const VarsAffine conditioned_rotation_;
  const VarsAffine rotation_;
  const double s_, mu1_;
  const double fopt_;
  std::vector<double> x_hat_, z1_, z2_; /**< @brief The transformed decision values */
};

/**
 * @brief The Lunacek bi-Rastrigin function f24 of the bbob-largescale suite, with permuted block rotations.
 *
 * f(x) = (10 (n - Σ cos(2π z_i)) + min(f_sphere(x^ - μ0), d n + s f_sphere(x^ - μ1))) / n + 10^4 f_pen(x) + fopt
 * with x^ = 2 sign x, the random signs sign, z = Q Λ^100 R (x^ - μ0) and the permuted block rotations R and Q.
 */
class BbobLunacekBiRastriginPermblockdiag final : public BbobProblem<BbobLunacekBiRastriginPermblockdiag> {
public:
  BbobLunacekBiRastriginPermblockdiag(const size_t function, const size_t dimension, const size_t instance,
                                      const long rseed, const char *problem_id_template,
                                      const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "5-weakly-structured", problem_id_template, problem_name_template),
        signs_(RandomSigns(rseed, dimension)), s_(1. - 0.5 / (sqrt(static_cast<double>(dimension + 20)) - 4.1)),
        mu1_(-sqrt((MU0 * MU0 - D) / s_)), shift_mu0_(std::vector<double>(dimension, MU0)),
        shift_mu1_(std::vector<double>(dimension, mu1_)),
        rotation1_(dimension, rseed + 3000000, rseed, rseed + 4000000), conditioning_(100.0, dimension),
        rotation2_(dimension, rseed + 5000000, rseed + 1000000, rseed + 6000000),
        fopt_(Bbob2009ComputeFopt(function, instance)), x_hat_(dimension), z1_(dimension), z2_(dimension) {
    /* The optimum is where x^ = μ0, map it back through the transformations */
    rotation2_.Inverse(best_parameter_, z1_);
    conditioning_.Inverse(z1_, z1_);
    rotation1_.Inverse(z1_, z2_);
    shift_mu0_.Inverse(z2_, z1_);
    for (size_t i = 0; i < dimension; ++i)
      best_parameter_[i] = z1_[i] / (2.0 * signs_[i]);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = BoundaryPenalty(x, smallest_values_of_interest_, largest_values_of_interest_);
    for (size_t i = 0; i < number_of_variables_; ++i)
      x_hat_[i] = 2.0 * signs_[i] * x[i];

    shift_mu0_.Apply(x_hat_, z1_);
    const double y0 = SphereRaw(z1_);
    rotation1_.Apply(z1_, z2_);
    conditioning_.Apply(z2_, z2_);
    rotation2_.Apply(z2_, z1_);
    shift_mu1_.Apply(x_hat_, z2_);
    const double y1 = SphereRaw(z2_) * s_ + D * static_cast<double>(number_of_variables_);

    double y = static_cast<double>(number_of_variables_);
    for (size_t i = 0; i < number_of_variables_; i++) {
      y -= cos(2 * coco_pi * z1_[i]);
    }
    y *= 10.0;
    y += DoubleMin(y0, y1);
    y *= NormByDimFactor(number_of_variables_);
    y += PENALTY_FACTOR * penalty;
    return y + fopt_;
  }

private:
  static constexpr double MU0 = 2.5;
  static constexpr double D = 1.0;
  static constexpr double PENALTY_FACTOR = 1e4;

  /** @brief Returns ±1 with random signs. */
  static std::vector<double> RandomSigns(const long rseed, const size_t dimension) {
    std::vector<double> signs(dimension);
    Bbob2009Gauss(signs, rseed);
    for (size_t i = 0; i < dimension; i++) {
      signs[i] = signs[i] < 0.0 ? -1.0 : 1.0;
    }
    return signs;
  }

  const std::vector<double> signs_;
  const double s_, mu1_;
  const VarsShift shift_mu0_; /**< @brief Subtracts μ0 */
  const VarsShift shift_mu1_; /**< @brief Subtracts μ1 */
  const VarsPermblockdiag rotation1_;
  const VarsConditioning conditioning_;
  const VarsPermblockdiag rotation2_;
  const double fopt_;
  std::vector<double> x_hat_, z1_, z2_; /**< @brief The transformed decision values */
};

} // namespace coco
