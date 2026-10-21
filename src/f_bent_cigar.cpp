/**
 * @file f_bent_cigar.cpp
 * @brief Implementation of the bent cigar function and problem.
 */

#include <assert.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_affine.cpp"
#include "transform_vars_shift.cpp"

namespace coco {

/**
 * @brief Implements the bent cigar function without connections to any COCO
 * structures.
 */
static double BentCigarRaw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  static const double condition = 1.0e6;

  double result = x[0] * x[0];
  for (size_t i = 1; i < number_of_variables; ++i) {
    result += condition * x[i] * x[i];
  }
  return result;
}

/**
 * @brief The bent cigar function.
 */
class BentCigar final : public Problem {
public:
  explicit BentCigar(const size_t number_of_variables) : Problem(number_of_variables, 1, 0) {
    InitFromScalars("bent cigar function", -5.0, 5.0, 0.0);
    SetId("%s_d%02zu", "bent_cigar", number_of_variables);

    /* Compute best solution */
    EvaluateFunction(best_parameter_, best_value_);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives_ == 1);
    y[0] = BentCigarRaw(x);
    assert(y[0] + 1e-13 >= best_value_[0]);
  }

  /**
   * @brief Evaluates the gradient of the bent cigar function.
   */
  void EvaluateGradient(std::span<const double> x, std::span<double> y) override {

    static const double condition = 1.0e6;

    y[0] = 2.0 * x[0];
    for (size_t i = 1; i < number_of_variables_; ++i)
      y[i] = 2.0 * condition * x[i];
  }
};

/**
 * @brief The BBOB bent cigar function f12, also in the bbob-boxed suite.
 *
 * f(x) = f_bent_cigar(R T_asy^0.5(R (x - xopt))) + fopt.
 */
class BbobBentCigar final : public BbobProblem<BbobBentCigar> {
public:
  BbobBentCigar(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                const BbobVariant variant, const char *problem_id_template, const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "3-ill-conditioned", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed + 1000000, dimension, variant)),
        rotation_(Bbob2009ComputeRotation(rseed + 1000000, dimension), std::vector<double>(dimension)),
        asymmetric_(0.5, dimension), fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension), z2_(dimension) {
    /* The optimum of the bent cigar function is the origin, map it back through the transformations */
    rotation_.Inverse(best_parameter_, z1_);
    asymmetric_.Inverse(z1_, z1_);
    rotation_.Inverse(z1_, z2_);
    shift_.Inverse(z2_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    shift_.Apply(x, z1_);
    rotation_.Apply(z1_, z2_);
    asymmetric_.Apply(z2_, z2_);
    rotation_.Apply(z2_, z1_);
    return BentCigarRaw(z1_) + fopt_;
  }

private:
  const VarsShift shift_;
  const VarsAffine rotation_;
  const VarsAsymmetric asymmetric_;
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

/**
 * @brief Creates the bent cigar problem for the constrained BBOB suite.
 */
static ProblemPtr BentCigarConsBbobProblemAllocate(const size_t function, const size_t dimension, const size_t instance,
                                                   const long rseed, const BbobVariant variant,
                                                   const char *problem_id_template, const char *problem_name_template) {

  ProblemPtr problem;

  Matrix M(dimension, dimension);
  std::vector<double> b(dimension);

  std::vector<double> xopt(dimension);
  double fopt = Bbob2009ComputeFopt(function, instance);
  BbobComputeXopt(xopt, rseed + 1000000, variant);

  const Matrix rot1 = Bbob2009ComputeRotation(rseed + 1000000, dimension);
  M = rot1;

  problem = std::make_unique<BentCigar>(dimension);
  problem = std::make_unique<TransformObjShift>(std::move(problem), fopt);
  problem = std::make_unique<TransformVarsAffine>(std::move(problem), M, b);
  problem = std::make_unique<TransformVarsShift>(std::move(problem), xopt, 0);

  problem->SetId(problem_id_template, function, instance, dimension);
  problem->SetName(problem_name_template, function, instance, dimension);
  problem->SetType("3-ill-conditioned");

  return problem;
}

} // namespace coco
