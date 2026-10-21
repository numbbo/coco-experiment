/**
 * @file f_ellipsoid.cpp
 * @brief Implementation of the ellipsoid function and problem.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_affine.cpp"

namespace coco {

/**
 * @brief The ellipsoid function.
 */
class Ellipsoid final : public Problem {
public:
  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double Raw(std::span<const double> x) const;

  Ellipsoid(const size_t number_of_variables, const double conditioning) : Problem(number_of_variables, 1, 0) {
    InitFromScalars("ellipsoid function", -5.0, 5.0, 0.0);
    SetId("%s_d%02zu", "ellipsoid", number_of_variables);

    scales_.resize(number_of_variables);
    for (size_t i = 0; i < number_of_variables; ++i) {
      const double exponent = 1.0 * static_cast<double>(i) / (static_cast<double>(number_of_variables) - 1.0);
      scales_[i] = pow(conditioning, exponent);
    }

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

  /**
   * @brief Evaluates the gradient of the ellipsoid function.
   */
  void EvaluateGradient(std::span<const double> x, std::span<double> y) override {
    for (size_t i = 0; i < number_of_variables_; ++i)
      y[i] = 2.0 * scales_[i] * x[i];
  }

private:
  std::vector<double> scales_; /**< @brief conditioning^(i / (n - 1)) */
};

/**
 * @brief Implements the ellipsoid function without connections to any COCO
 * structures.
 */
double Ellipsoid::Raw(std::span<const double> x) const {

  double result = x[0] * x[0];
  for (size_t i = 1; i < number_of_variables_; ++i)
    result += scales_[i] * x[i] * x[i];

  return result;
}

/**
 * @brief The BBOB ellipsoid function f2, also in the bbob-boxed and bbob-largescale suites.
 *
 * f(x) = f_ellipsoid^1e6(T_osz(x - xopt)) + fopt, normalized by the dimension before fopt is added in the
 * large-scale suite.
 */
class BbobEllipsoid final : public BbobProblem<BbobEllipsoid> {
public:
  BbobEllipsoid(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                const BbobVariant variant, const char *problem_id_template, const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "1-separable", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)), ellipsoid_(dimension, 1.0e6),
        fopt_(Bbob2009ComputeFopt(function, instance)), normalize_(variant == BbobVariant::LARGESCALE), z_(dimension) {
    /* The optimum of the ellipsoid is the origin, map it back through the transformations */
    oscillate_.Inverse(best_parameter_, z_);
    shift_.Inverse(z_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    shift_.Apply(x, z_);
    oscillate_.Apply(z_, z_);
    double y = ellipsoid_.Raw(z_);
    if (normalize_)
      y *= NormByDimFactor(number_of_variables_);
    return y + fopt_;
  }

private:
  const VarsShift shift_;
  const VarsOscillate oscillate_{};
  const Ellipsoid ellipsoid_;
  const double fopt_;
  const bool normalize_;
  std::vector<double> z_; /**< @brief The transformed decision values */
};

/**
 * @brief The BBOB rotated ellipsoid function f10, also in the bbob-boxed and bbob-noisy suites.
 *
 * f(x) = f_ellipsoid^conditioning(T_osz(R (x - xopt))) + fopt.
 */
class BbobEllipsoidRotated final : public BbobProblem<BbobEllipsoidRotated> {
public:
  BbobEllipsoidRotated(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                       const BbobVariant variant, const double conditioning, const char *problem_id_template,
                       const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "3-ill-conditioned", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)),
        rotation_(Bbob2009ComputeRotation(rseed + 1000000, dimension), std::vector<double>(dimension)),
        ellipsoid_(dimension, conditioning), fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension),
        z2_(dimension) {
    /* The optimum of the ellipsoid is the origin, map it back through the transformations */
    oscillate_.Inverse(best_parameter_, z1_);
    rotation_.Inverse(z1_, z2_);
    shift_.Inverse(z2_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    shift_.Apply(x, z1_);
    rotation_.Apply(z1_, z2_);
    oscillate_.Apply(z2_, z2_);
    return ellipsoid_.Raw(z2_) + fopt_;
  }

private:
  const VarsShift shift_;
  const VarsAffine rotation_;
  const VarsOscillate oscillate_{};
  const Ellipsoid ellipsoid_;
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

/**
 * @brief The rotated ellipsoid function f10 of the bbob-largescale suite, with a permuted block rotation.
 *
 * f(x) = f_ellipsoid^1e6(T_osz(R (x - xopt))) / n + fopt with the permuted block rotation R.
 */
class BbobEllipsoidPermblockdiag final : public BbobProblem<BbobEllipsoidPermblockdiag> {
public:
  BbobEllipsoidPermblockdiag(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                             const BbobVariant variant, const char *problem_id_template,
                             const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "3-ill-conditioned", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)),
        rotation_(dimension, rseed + 2000000, rseed + 1000000, rseed + 3000000), ellipsoid_(dimension, 1.0e6),
        fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension), z2_(dimension) {
    /* The optimum of the ellipsoid is the origin, map it back through the transformations */
    oscillate_.Inverse(best_parameter_, z1_);
    rotation_.Inverse(z1_, z2_);
    shift_.Inverse(z2_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    shift_.Apply(x, z1_);
    rotation_.Apply(z1_, z2_);
    oscillate_.Apply(z2_, z2_);
    return ellipsoid_.Raw(z2_) * NormByDimFactor(number_of_variables_) + fopt_;
  }

private:
  const VarsShift shift_;
  const VarsPermblockdiag rotation_;
  const VarsOscillate oscillate_{};
  const Ellipsoid ellipsoid_;
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

/**
 * @brief Creates the ellipsoid problem for the constrained BBOB suite
 */
static ProblemPtr EllipsoidConsBbobProblemAllocate(const size_t function, const size_t dimension, const size_t instance,
                                                   const long rseed, const BbobVariant variant,
                                                   const char *problem_id_template, const char *problem_name_template) {
  ProblemPtr problem;

  std::vector<double> xopt(dimension);
  double fopt = Bbob2009ComputeFopt(function, instance);
  BbobComputeXopt(xopt, rseed, variant);

  problem = std::make_unique<Ellipsoid>(dimension, 1.0e6);

  /* TODO (NH): fopt -= problem->evaluate(all_zeros(dimension)) */
  problem = std::make_unique<TransformVarsShift>(std::move(problem), xopt, 0);
  problem = std::make_unique<TransformObjShift>(std::move(problem), fopt);

  problem->SetId(problem_id_template, function, instance, dimension);
  problem->SetName(problem_name_template, function, instance, dimension);
  problem->SetType("1-separable");

  return problem;
}

/**
 * @brief Creates the rotated ellipsoid problem for the constrained
 *        BBOB suite
 */
static ProblemPtr EllipsoidRotatedConsBbobProblemAllocate(const size_t function, const size_t dimension,
                                                          const size_t instance, const long rseed,
                                                          const BbobVariant variant, const char *problem_id_template,
                                                          const char *problem_name_template) {
  ProblemPtr problem;

  Matrix M(dimension, dimension);
  std::vector<double> b(dimension);

  std::vector<double> xopt(dimension);
  BbobComputeXopt(xopt, rseed, variant);
  double fopt = Bbob2009ComputeFopt(function, instance);

  const Matrix rot1 = Bbob2009ComputeRotation(rseed + 1000000, dimension);
  M = rot1;

  problem = std::make_unique<Ellipsoid>(dimension, 1.0e6);
  problem = std::make_unique<TransformVarsAffine>(std::move(problem), M, b);
  problem = std::make_unique<TransformVarsShift>(std::move(problem), xopt, 0);
  problem = std::make_unique<TransformObjShift>(std::move(problem), fopt);

  problem->SetId(problem_id_template, function, instance, dimension);
  problem->SetName(problem_name_template, function, instance, dimension);
  problem->SetType("3-ill-conditioned");

  return problem;
}

} // namespace coco
