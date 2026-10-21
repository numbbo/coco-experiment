/**
 * @file f_linear_slope.cpp
 * @brief Implementation of the linear slope function and problem.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_shift.cpp"

namespace coco {

/**
 * @brief Implements the linear slope function without connections to any COCO
 * structures.
 */
static double LinearSlopeRaw(std::span<const double> x, std::span<const double> best_parameter,
                             std::span<const double> slopes) {
  double result = 0.0;

  assert(best_parameter.size() == x.size() && slopes.size() == x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    const double si = best_parameter[i] > 0.0 ? slopes[i] : -slopes[i];
    /* boundary handling */
    if (x[i] * best_parameter[i] < 25.0) {
      result += 5.0 * fabs(si) - si * x[i];
    } else {
      result += 5.0 * fabs(si) - si * best_parameter[i];
    }
  }

  return result;
}

/**
 * @brief The linear slope function.
 */
class LinearSlope final : public Problem {
public:
  LinearSlope(const size_t number_of_variables, std::span<const double> best_parameter)
      : Problem(number_of_variables, 1, 0), slopes_(number_of_variables) {
    assert(best_parameter.size() == number_of_variables);
    static const double alpha = 100.0;
    for (size_t i = 0; i < number_of_variables; ++i)
      slopes_[i] = pow(sqrt(alpha), static_cast<double>(i) / (static_cast<double>(number_of_variables) - 1));

    /* best_parameter will be overwritten below */
    InitFromScalars("linear slope function", -5.0, 5.0, 0.0);
    SetId("%s_d%02zu", "linear_slope", number_of_variables);

    /* Compute best solution */
    for (size_t i = 0; i < number_of_variables; ++i) {
      if (best_parameter[i] < 0.0) {
        this->best_parameter_[i] = smallest_values_of_interest_[i];
      } else {
        this->best_parameter_[i] = largest_values_of_interest_[i];
      }
    }
    EvaluateFunction(this->best_parameter_, best_value_);
  }

  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double Raw(std::span<const double> x) const { return LinearSlopeRaw(x, best_parameter_, slopes_); }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives_ == 1);
    y[0] = Raw(x);
    assert(y[0] + 1e-13 >= best_value_[0]);
  }

  /**
   * @brief Evaluates the gradient of the linear slope function.
   */
  void EvaluateGradient(std::span<const double>, std::span<double> y) override {
    for (size_t i = 0; i < number_of_variables_; ++i)
      y[i] = best_parameter_[i] > 0.0 ? -slopes_[i] : slopes_[i];
  }

private:
  std::vector<double> slopes_; /**< @brief |s_i| = sqrt(alpha)^(i / (n - 1)) */
};

/**
 * @brief The BBOB linear slope function f5, also in the bbob-boxed and bbob-largescale suites.
 *
 * f(x) = f_linear_slope(x) + fopt with the optimum at the corner of the region of interest given by the signs of
 * xopt, normalized by the dimension before fopt is added in the large-scale suite.
 */
class BbobLinearSlope final : public BbobProblem<BbobLinearSlope> {
public:
  BbobLinearSlope(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                  const BbobVariant variant, const char *problem_id_template, const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "1-separable", problem_id_template, problem_name_template),
        linear_slope_(dimension, BbobXopt(rseed, dimension, variant)), fopt_(Bbob2009ComputeFopt(function, instance)),
        normalize_(variant == BbobVariant::LARGESCALE) {
    best_parameter_ = linear_slope_.best_parameter();
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    double y = linear_slope_.Raw(x);
    if (normalize_)
      y *= NormByDimFactor(number_of_variables_);
    return y + fopt_;
  }

private:
  const LinearSlope linear_slope_;
  const double fopt_;
  const bool normalize_;
};

/**
 * @brief Creates the linear slope problem for the constrained BBOB suite.
 */
static ProblemPtr LinearSlopeConsBbobProblemAllocate(const size_t function, const size_t dimension,
                                                     const size_t instance, const long rseed,
                                                     const char *problem_id_template,
                                                     const char *problem_name_template) {
  ProblemPtr problem;

  std::vector<double> xopt(dimension);
  BbobComputeXopt(xopt, rseed, BbobVariant::BBOB);
  double fopt = Bbob2009ComputeFopt(function, instance);

  problem = std::make_unique<LinearSlope>(dimension, xopt);
  problem = std::make_unique<TransformObjShift>(std::move(problem), fopt);

  problem->SetId(problem_id_template, function, instance, dimension);
  problem->SetName(problem_name_template, function, instance, dimension);
  problem->SetType("1-separable");

  return problem;
}

} // namespace coco
