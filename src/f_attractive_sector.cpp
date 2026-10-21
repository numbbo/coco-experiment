/**
 * @file f_attractive_sector.cpp
 * @brief Implementation of the attractive sector function and problem.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

namespace coco {

/**
 * @brief The attractive sector function.
 */
class AttractiveSector final : public Problem {
public:
  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double Raw(std::span<const double> x) const;

  AttractiveSector(const size_t number_of_variables, std::span<const double> xopt)
      : Problem(number_of_variables, 1, 0) {
    assert(xopt.size() == number_of_variables);
    InitFromScalars("attractive sector function", -5.0, 5.0, 0.0);
    SetId("%s_d%02zu", "attractive_sector", number_of_variables);

    this->xopt_.assign(xopt.begin(), xopt.end());

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
  std::vector<double> xopt_;
};

/**
 * @brief Implements the attractive sector function without connections to any
 * COCO structures.
 */
double AttractiveSector::Raw(std::span<const double> x) const {

  double result = 0.0;
  for (size_t i = 0; i < number_of_variables_; ++i) {
    if (xopt_[i] * x[i] > 0.0) {
      result += 100.0 * 100.0 * x[i] * x[i];
    } else {
      result += x[i] * x[i];
    }
  }
  return result;
}

/**
 * @brief The BBOB attractive sector function f6, also in the bbob-boxed suite.
 *
 * f(x) = T_osz(f_attractive_sector(Q Λ^10 R (x - xopt)))^0.9 + fopt.
 */
class BbobAttractiveSector final : public BbobProblem<BbobAttractiveSector> {
public:
  BbobAttractiveSector(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                       const BbobVariant variant, const char *problem_id_template, const char *problem_name_template)
      : BbobAttractiveSector(function, dimension, instance, rseed, BbobXopt(rseed, dimension, variant),
                             problem_id_template, problem_name_template) {}

  double operator()(std::span<const double> x) {
    shift_.Apply(x, z1_);
    conditioned_rotation_.Apply(z1_, z2_);
    const double y = pow(OscillateObjective(sector_.Raw(z2_)), 0.9);
    return y + fopt_;
  }

private:
  BbobAttractiveSector(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                       const std::vector<double> &xopt, const char *problem_id_template,
                       const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "2-moderate", problem_id_template, problem_name_template),
        shift_(xopt),
        conditioned_rotation_(ConditionedRotationMatrix(rseed, dimension), std::vector<double>(dimension)),
        sector_(dimension, xopt), fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension), z2_(dimension) {
    /* The optimum of the attractive sector function is the origin, map it back through the transformations */
    conditioned_rotation_.Inverse(best_parameter_, z1_);
    shift_.Inverse(z1_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  /** @brief Returns Q Λ^10 R with the rotations Q and R. */
  static Matrix ConditionedRotationMatrix(const long rseed, const size_t dimension) {
    const Matrix rot1 = Bbob2009ComputeRotation(rseed + 1000000, dimension);
    const Matrix rot2 = Bbob2009ComputeRotation(rseed, dimension);
    std::vector<double> scales(dimension);
    for (size_t k = 0; k < dimension; ++k) {
      const double exponent = 1.0 * static_cast<int>(k) / (static_cast<double>(dimension) - 1.0);
      scales[k] = pow(sqrt(10.0), exponent);
    }
    return MatrixProductWithDiagonal(rot1, scales, rot2);
  }

  const VarsShift shift_;
  const VarsAffine conditioned_rotation_;
  const AttractiveSector sector_;
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

/**
 * @brief The attractive sector function f6 of the bbob-largescale suite, with permuted block rotations.
 *
 * f(x) = T_osz(f_attractive_sector(Q Λ^10 R (x - xopt)) / n)^0.9 + fopt with the permuted block rotations R and Q.
 */
class BbobAttractiveSectorPermblockdiag final : public BbobProblem<BbobAttractiveSectorPermblockdiag> {
public:
  BbobAttractiveSectorPermblockdiag(const size_t function, const size_t dimension, const size_t instance,
                                    const long rseed, const BbobVariant variant, const char *problem_id_template,
                                    const char *problem_name_template)
      : BbobAttractiveSectorPermblockdiag(function, dimension, instance, rseed, BbobXopt(rseed, dimension, variant),
                                          problem_id_template, problem_name_template) {}

  double operator()(std::span<const double> x) {
    shift_.Apply(x, z1_);
    rotation1_.Apply(z1_, z2_);
    conditioning_.Apply(z2_, z2_);
    rotation2_.Apply(z2_, z1_);
    double y = sector_.Raw(z1_) * NormByDimFactor(number_of_variables_);
    y = pow(OscillateObjective(y), 0.9);
    return y + fopt_;
  }

private:
  BbobAttractiveSectorPermblockdiag(const size_t function, const size_t dimension, const size_t instance,
                                    const long rseed, const std::vector<double> &xopt, const char *problem_id_template,
                                    const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "2-moderate", problem_id_template, problem_name_template),
        shift_(xopt), rotation1_(dimension, rseed + 2000000, rseed, rseed + 3000000), conditioning_(10.0, dimension),
        rotation2_(dimension, rseed + 4000000, rseed + 1000000, rseed + 5000000), sector_(dimension, xopt),
        fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension), z2_(dimension) {
    /* The optimum of the attractive sector function is the origin, map it back through the transformations */
    rotation2_.Inverse(best_parameter_, z1_);
    conditioning_.Inverse(z1_, z1_);
    rotation1_.Inverse(z1_, z2_);
    shift_.Inverse(z2_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  const VarsShift shift_;
  const VarsPermblockdiag rotation1_;
  const VarsConditioning conditioning_;
  const VarsPermblockdiag rotation2_;
  const AttractiveSector sector_;
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

} // namespace coco
