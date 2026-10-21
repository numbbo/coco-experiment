/**
 * @file f_discus_generalized.cpp
 * @brief Implementation of the discus function and problem.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

namespace coco {

/**
 * @brief The generalized discus function.
 */
class DiscusGeneralized final : public Problem {
public:
  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double Raw(std::span<const double> x) const;

  DiscusGeneralized(const size_t number_of_variables, size_t proportion_short_axes_denom)
      : Problem(number_of_variables, 1, 0) {
    InitFromScalars("generalized discus function", -5.0, 5.0, 0.0);
    SetId("%s_d%04zu", "discus_generalized", number_of_variables);
    this->proportion_short_axes_denom_ = proportion_short_axes_denom;

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
  size_t proportion_short_axes_denom_{};
};

/**
 * @brief Implements the generalized discus function without connections to any
 * COCO structures.
 */
double DiscusGeneralized::Raw(std::span<const double> x) const {

  static const double condition = 1.0e6;
  double result = 0;
  size_t nb_short_axes = number_of_variables_ / proportion_short_axes_denom_;
  if (number_of_variables_ % proportion_short_axes_denom_ != 0) {
    nb_short_axes += 1;
  }

  for (size_t i = 0; i < nb_short_axes; ++i) {
    result += x[i] * x[i];
  }
  result *= condition;

  for (size_t i = nb_short_axes; i < number_of_variables_; ++i) {
    result += x[i] * x[i];
  }

  return result;
}

/**
 * @brief The discus function f11 of the bbob-largescale suite, generalized to n / 40 short axes and with a
 * permuted block rotation.
 *
 * f(x) = f_discus_generalized(T_osz(R (x - xopt))) / n + fopt with the permuted block rotation R.
 */
class BbobDiscusGeneralizedPermblockdiag final : public BbobProblem<BbobDiscusGeneralizedPermblockdiag> {
public:
  BbobDiscusGeneralizedPermblockdiag(const size_t function, const size_t dimension, const size_t instance,
                                     const long rseed, const BbobVariant variant, const char *problem_id_template,
                                     const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "3-ill-conditioned", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)),
        rotation_(dimension, rseed + 2000000, rseed + 1000000, rseed + 3000000), discus_(dimension, 40),
        fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension), z2_(dimension) {
    /* The optimum of the discus function is the origin, map it back through the transformations */
    oscillate_.Inverse(best_parameter_, z1_);
    rotation_.Inverse(z1_, z2_);
    shift_.Inverse(z2_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    shift_.Apply(x, z1_);
    rotation_.Apply(z1_, z2_);
    oscillate_.Apply(z2_, z2_);
    return discus_.Raw(z2_) * NormByDimFactor(number_of_variables_) + fopt_;
  }

private:
  const VarsShift shift_;
  const VarsPermblockdiag rotation_;
  const VarsOscillate oscillate_{};
  const DiscusGeneralized discus_;
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

} // namespace coco
