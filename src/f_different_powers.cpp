/**
 * @file f_different_powers.cpp
 * @brief Implementation of the different powers function and problem.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_affine.cpp"
#include "transform_vars_shift.cpp"

namespace coco {

/**
 * @brief Implements the different powers function without connections to any
 * COCO structures.
 */
static double DifferentPowersRaw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  double sum = 0.0;

  for (size_t i = 0; i < number_of_variables; ++i) {
    double exponent = 2.0 + (4.0 * static_cast<double>(i)) / (static_cast<double>(number_of_variables) - 1.0);
    sum += pow(fabs(x[i]), exponent);
  }
  double result = sqrt(sum);

  return result;
}

/**
 * @brief Implements the sign function.
 */
double Sign(double x) {

  if (x > 0)
    return 1;
  if (x < 0)
    return -1;
  return 0;
}

/**
 * @brief The different powers function.
 */
class DifferentPowers final : public Problem {
public:
  explicit DifferentPowers(const size_t number_of_variables) : Problem(number_of_variables, 1, 0) {
    InitFromScalars("different powers function", -5.0, 5.0, 0.0);
    SetId("%s_d%02zu", "different_powers", number_of_variables);

    /* Compute best solution */
    EvaluateFunction(best_parameter_, best_value_);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives_ == 1);
    y[0] = DifferentPowersRaw(x);
    assert(y[0] + 1e-13 >= best_value_[0]);
  }

  /**
   * @brief Evaluates the gradient of the function "different powers".
   */
  void EvaluateGradient(std::span<const double> x, std::span<double> y) override {

    double sum = 0.0;
    double aux;

    for (size_t i = 0; i < number_of_variables_; ++i) {
      aux = 2.0 + (4.0 * static_cast<double>(i)) / (static_cast<double>(number_of_variables_) - 1.0);
      sum += pow(fabs(x[i]), aux);
    }

    for (size_t i = 0; i < number_of_variables_; ++i) {
      aux = 2.0 + (4.0 * static_cast<double>(i)) / (static_cast<double>(number_of_variables_) - 1.0);
      y[i] = 0.5 * (aux) / (sum);
      aux -= 1.0;
      y[i] *= pow(fabs(x[i]), aux) * Sign(x[i]);
    }
  }
};

/**
 * @brief The BBOB different powers function f14, also in the bbob-boxed and bbob-noisy suites.
 *
 * f(x) = f_different_powers(R (x - xopt)) + fopt.
 */
class BbobDifferentPowers final : public BbobProblem<BbobDifferentPowers> {
public:
  BbobDifferentPowers(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                      const BbobVariant variant, const char *problem_id_template, const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "3-ill-conditioned", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)),
        rotation_(Bbob2009ComputeRotation(rseed + 1000000, dimension), std::vector<double>(dimension)),
        fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension), z2_(dimension) {
    /* The optimum of the different powers function is the origin, map it back through the transformations */
    rotation_.Inverse(best_parameter_, z1_);
    shift_.Inverse(z1_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    shift_.Apply(x, z1_);
    rotation_.Apply(z1_, z2_);
    return DifferentPowersRaw(z2_) + fopt_;
  }

private:
  const VarsShift shift_;
  const VarsAffine rotation_;
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

/**
 * @brief Creates the bbob-constrained different powers problem.
 */
static ProblemPtr DifferentPowersBbobConstrainedProblemAllocate(const size_t function, const size_t dimension,
                                                                const size_t instance, const long rseed,
                                                                const BbobVariant variant,
                                                                const char *problem_id_template,
                                                                const char *problem_name_template) {
  /* Different powers function used in bbob-constrained test suite.
   * In this version, the (unconstrained) optimum, xopt, is set to
   * a distance of 1e-2 to the origin. By doing so, the optimum of
   * the constrained problem is at a "reasonable" distance from
   * the unconstrained one and, hence, the constrained problem is not too easy.
   */

  double fopt, result;
  ProblemPtr problem;

  Matrix M(dimension, dimension);
  std::vector<double> b(dimension);

  std::vector<double> xopt(dimension);
  fopt = Bbob2009ComputeFopt(function, instance);
  BbobComputeXopt(xopt, rseed, variant);

  /* Compute Euclidean norm of xopt */
  /* CAVEAT: this implementation is not ideal for large dimensions */
  result = 0.0;
  for (size_t i = 0; i < dimension; ++i) {
    result += xopt[i] * xopt[i];
  }
  result = sqrt(result);

  /* Scale xopt such that the distance to the origin is 1e-2 */
  for (size_t i = 0; i < dimension; ++i) {
    xopt[i] *= 1e-2 / result;
  }

  const Matrix rot1 = Bbob2009ComputeRotation(rseed + 1000000, dimension);
  M = rot1;

  problem = std::make_unique<DifferentPowers>(dimension);
  problem = std::make_unique<TransformObjShift>(std::move(problem), fopt);
  problem = std::make_unique<TransformVarsAffine>(std::move(problem), M, b);
  problem = std::make_unique<TransformVarsShift>(std::move(problem), xopt, 0);

  problem->SetId(problem_id_template, function, instance, dimension);
  problem->SetName(problem_name_template, function, instance, dimension);
  problem->SetType("3-ill-conditioned");

  return problem;
}

/**
 * @brief The different powers function f14 of the bbob-largescale suite, with a permuted block rotation.
 *
 * f(x) = f_different_powers(R (x - xopt)) / n + fopt with the permuted block rotation R.
 */
class BbobDifferentPowersPermblockdiag final : public BbobProblem<BbobDifferentPowersPermblockdiag> {
public:
  BbobDifferentPowersPermblockdiag(const size_t function, const size_t dimension, const size_t instance,
                                   const long rseed, const BbobVariant variant, const char *problem_id_template,
                                   const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "3-ill-conditioned", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)),
        rotation_(dimension, rseed + 2000000, rseed + 1000000, rseed + 3000000),
        fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension), z2_(dimension) {
    /* The optimum of the different powers function is the origin, map it back through the transformations */
    rotation_.Inverse(best_parameter_, z1_);
    shift_.Inverse(z1_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    shift_.Apply(x, z1_);
    rotation_.Apply(z1_, z2_);
    return DifferentPowersRaw(z2_) * NormByDimFactor(number_of_variables_) + fopt_;
  }

private:
  const VarsShift shift_;
  const VarsPermblockdiag rotation_;
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

} // namespace coco
