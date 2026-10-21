/**
 * @file f_rastrigin.cpp
 * @brief Implementation of the Rastrigin function and problem.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "coco_utilities.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_affine.cpp"
#include "transform_vars_shift.cpp"
#include "bbob_problem.cpp"
#include "transform_functions.cpp"

namespace coco {

/**
 * @brief Implements the Rastrigin function without connections to any COCO
 * structures.
 */
static double RastriginRaw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  size_t i = 0;
  double sum1 = 0.0, sum2 = 0.0;

  for (i = 0; i < number_of_variables; ++i) {
    sum1 += cos(coco_two_pi * x[i]);
    sum2 += x[i] * x[i];
  }
  if (IsInf(sum2)) /* cos(inf) -> nan */
    return sum2;
  double result = 10.0 * (static_cast<double>(number_of_variables) - sum1) + sum2;

  return result;
}

/**
 * @brief The Rastrigin function.
 */
class Rastrigin final : public Problem {
public:
  explicit Rastrigin(const size_t number_of_variables) : Problem(number_of_variables, 1, 0) {
    InitFromScalars("Rastrigin function", -5.0, 5.0, 0.0);
    /* TODO: make sure the gradient is computed correctly for the rotated Rastrigin */
    SetId("%s_d%02zu", "rastrigin", number_of_variables);

    /* Compute best solution */
    EvaluateFunction(best_parameter_, best_value_);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives_ == 1);
    y[0] = RastriginRaw(x);
    assert(y[0] + 1e-13 >= best_value_[0]);
  }

  /**
   * @brief Evaluates the gradient of the raw Rastrigin function.
   */
  void EvaluateGradient(std::span<const double> x, std::span<double> y) override {

    for (size_t i = 0; i < number_of_variables_; ++i) {
      y[i] = 2.0 * (10. * coco_pi * sin(coco_two_pi * x[i]) + x[i]);
    }
  }
};

/**
 * @brief The BBOB Rastrigin function f3, also in the bbob-boxed and bbob-largescale suites.
 *
 * f(x) = f_rastrigin(Λ^10 T_asy^0.2(T_osz(x - xopt))) + fopt, normalized by the dimension before fopt is added in
 * the large-scale suite.
 */
class BbobRastrigin final : public BbobProblem<BbobRastrigin> {
public:
  BbobRastrigin(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                const BbobVariant variant, const char *problem_id_template, const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "1-separable", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)), asymmetric_(0.2, dimension), conditioning_(10.0, dimension),
        fopt_(Bbob2009ComputeFopt(function, instance)), normalize_(variant == BbobVariant::LARGESCALE), z_(dimension) {
    /* The optimum of the Rastrigin function is the origin, map it back through the transformations */
    conditioning_.Inverse(best_parameter_, z_);
    asymmetric_.Inverse(z_, z_);
    oscillate_.Inverse(z_, z_);
    shift_.Inverse(z_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    shift_.Apply(x, z_);
    oscillate_.Apply(z_, z_);
    asymmetric_.Apply(z_, z_);
    conditioning_.Apply(z_, z_);
    double y = RastriginRaw(z_);
    if (normalize_)
      y *= NormByDimFactor(number_of_variables_);
    return y + fopt_;
  }

private:
  const VarsShift shift_;
  const VarsOscillate oscillate_{};
  const VarsAsymmetric asymmetric_;
  const VarsConditioning conditioning_;
  const double fopt_;
  const bool normalize_;
  std::vector<double> z_; /**< @brief The transformed decision values */
};

/**
 * @brief The BBOB rotated Rastrigin function f15, also in the bbob-boxed suite.
 *
 * f(x) = f_rastrigin(R Λ^10 Q T_asy^0.2(T_osz(R (x - xopt)))) + fopt.
 */
class BbobRastriginRotated final : public BbobProblem<BbobRastriginRotated> {
public:
  BbobRastriginRotated(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                       const BbobVariant variant, const char *problem_id_template, const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "4-multi-modal", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)),
        rotation_(Bbob2009ComputeRotation(rseed + 1000000, dimension), std::vector<double>(dimension)),
        asymmetric_(0.2, dimension),
        conditioned_rotation_(ConditionedRotationMatrix(rseed, dimension), std::vector<double>(dimension)),
        fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension), z2_(dimension) {
    /* The optimum of the Rastrigin function is the origin, map it back through the transformations */
    conditioned_rotation_.Inverse(best_parameter_, z1_);
    asymmetric_.Inverse(z1_, z1_);
    oscillate_.Inverse(z1_, z1_);
    rotation_.Inverse(z1_, z2_);
    shift_.Inverse(z2_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    shift_.Apply(x, z1_);
    rotation_.Apply(z1_, z2_);
    oscillate_.Apply(z2_, z2_);
    asymmetric_.Apply(z2_, z2_);
    conditioned_rotation_.Apply(z2_, z1_);
    return RastriginRaw(z1_) + fopt_;
  }

private:
  /** @brief Returns R Λ^10 Q with the rotations R and Q. */
  static Matrix ConditionedRotationMatrix(const long rseed, const size_t dimension) {
    const Matrix rot1 = Bbob2009ComputeRotation(rseed + 1000000, dimension);
    const Matrix rot2 = Bbob2009ComputeRotation(rseed, dimension);
    std::vector<double> scales(dimension);
    for (size_t k = 0; k < dimension; ++k) {
      const double exponent = 1.0 * static_cast<int>(k) / (static_cast<double>(dimension) - 1.0);
      scales[k] = pow(sqrt(10), exponent);
    }
    return MatrixProductWithDiagonal(rot1, scales, rot2);
  }

  const VarsShift shift_;
  const VarsAffine rotation_;
  const VarsOscillate oscillate_{};
  const VarsAsymmetric asymmetric_;
  const VarsAffine conditioned_rotation_;
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

/**
 * @brief The rotated Rastrigin function f15 of the bbob-largescale suite, with permuted block rotations.
 *
 * f(x) = f_rastrigin(R Λ^10 Q T_asy^0.2(T_osz(R (x - xopt)))) / n + fopt with the permuted block rotations R and Q.
 */
class BbobRastriginRotatedPermblockdiag final : public BbobProblem<BbobRastriginRotatedPermblockdiag> {
public:
  BbobRastriginRotatedPermblockdiag(const size_t function, const size_t dimension, const size_t instance,
                                    const long rseed, const BbobVariant variant, const char *problem_id_template,
                                    const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "4-multi-modal", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)),
        rotation1_(dimension, rseed + 3000000, rseed + 1000000, rseed + 4000000), asymmetric_(0.2, dimension),
        rotation2_(dimension, rseed + 5000000, rseed, rseed + 6000000), conditioning_(10.0, dimension),
        fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension), z2_(dimension) {
    /* The optimum of the Rastrigin function is the origin, map it back through the transformations */
    rotation1_.Inverse(best_parameter_, z1_);
    conditioning_.Inverse(z1_, z1_);
    rotation2_.Inverse(z1_, z2_);
    asymmetric_.Inverse(z2_, z2_);
    oscillate_.Inverse(z2_, z2_);
    rotation1_.Inverse(z2_, z1_);
    shift_.Inverse(z1_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    shift_.Apply(x, z1_);
    rotation1_.Apply(z1_, z2_);
    oscillate_.Apply(z2_, z2_);
    asymmetric_.Apply(z2_, z2_);
    rotation2_.Apply(z2_, z1_);
    conditioning_.Apply(z1_, z1_);
    rotation1_.Apply(z1_, z2_);
    const double y = RastriginRaw(z2_);
    return y * NormByDimFactor(number_of_variables_) + fopt_;
  }

private:
  const VarsShift shift_;
  const VarsPermblockdiag rotation1_;
  const VarsOscillate oscillate_{};
  const VarsAsymmetric asymmetric_;
  const VarsPermblockdiag rotation2_;
  const VarsConditioning conditioning_;
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

/**
 * @brief Computes xopt for constrained Rastrigin (alternative to
 * Bbob2009ComputeXopt()) xopt is a vector of dim uniform random integers
 */
static void RastriginConsComputeXopt(std::span<double> xopt, const long rseed) {
  Bbob2009Unif(xopt, rseed);

  for (size_t i = 0; i < xopt.size(); ++i) {
    xopt[i] = 10 * xopt[i] - 5;
    xopt[i] = static_cast<int>(xopt[i]);
  }

  /* In case (0, ..., 0) is sampled, set xopt to a different value */
  if (VectorIsZero(xopt))
    for (size_t i = 0; i < xopt.size(); ++i) {
      xopt[i] = static_cast<int>(i % 9) - 4;
    }
}

/**
 * @brief Creates the Rastrigin problem for the constrained BBOB suite.
 */
static ProblemPtr RastriginConsBbobProblemAllocate(const size_t function, const size_t dimension, const size_t instance,
                                                   const long rseed, const char *problem_id_template,
                                                   const char *problem_name_template) {

  ProblemPtr problem;

  std::vector<double> xopt(dimension);
  RastriginConsComputeXopt(xopt, rseed);
  double fopt = Bbob2009ComputeFopt(function, instance);

  problem = std::make_unique<Rastrigin>(dimension);
  problem = std::make_unique<TransformVarsShift>(std::move(problem), xopt, 0);
  problem = std::make_unique<TransformObjShift>(std::move(problem), fopt);

  problem->SetId(problem_id_template, function, instance, dimension);
  problem->SetName(problem_name_template, function, instance, dimension);
  problem->SetType("1-separable");

  return problem;
}

/**
 * @brief Creates the rotated Rastrigin problem for the constrained BBOB suite.
 */
static ProblemPtr RastriginRotatedConsBbobProblemAllocate(const size_t function, const size_t dimension,
                                                          const size_t instance, const long rseed,
                                                          const char *problem_id_template,
                                                          const char *problem_name_template) {

  ProblemPtr problem;

  Matrix M(dimension, dimension);
  std::vector<double> b(dimension);

  const Matrix rot1 = Bbob2009ComputeRotation(rseed + 1000000, dimension);
  M = rot1;

  std::vector<double> xopt(dimension);
  RastriginConsComputeXopt(xopt, rseed);
  double fopt = Bbob2009ComputeFopt(function, instance);

  problem = std::make_unique<Rastrigin>(dimension);
  problem = std::make_unique<TransformVarsShift>(std::move(problem), xopt, 0);
  problem = std::make_unique<TransformVarsAffine>(std::move(problem), M,
                                                  b); /* Rotate after shifting so R dot xopt is also rotated */
  problem = std::make_unique<TransformObjShift>(std::move(problem), fopt);

  problem->SetId(problem_id_template, function, instance, dimension);
  problem->SetName(problem_name_template, function, instance, dimension);
  problem->SetType("1-separable");

  return problem;
}

} // namespace coco
