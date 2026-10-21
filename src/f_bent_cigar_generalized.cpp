/**
 * @file f_bent_cigar_generalized.cpp
 * @brief Implementation of the generalized bent cigar function and problem.
 */

#include <assert.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

namespace coco {

/**
 * @brief The generalized bent cigar function.
 */
class BentCigarGeneralized final : public Problem {
public:
  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double Raw(std::span<const double> x) const;

  BentCigarGeneralized(const size_t number_of_variables, size_t proportion_long_axes_denom)
      : Problem(number_of_variables, 1, 0) {
    InitFromScalars("generalized bent cigar function", -5.0, 5.0, 0.0);
    SetId("%s_d%04zu", "bent_cigar", number_of_variables);
    this->proportion_long_axes_denom_ = proportion_long_axes_denom;

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
  size_t proportion_long_axes_denom_{};
};

/**
 * @brief Implements the generalized bent cigar function without connections to
 * any COCO structures.
 */
double BentCigarGeneralized::Raw(std::span<const double> x) const {

  static const double condition = 1.0e6;
  double result = 0;
  size_t nb_long_axes = number_of_variables_ / proportion_long_axes_denom_;
  if (number_of_variables_ % proportion_long_axes_denom_ != 0) {
    nb_long_axes += 1;
  }
  for (size_t i = 0; i < nb_long_axes; ++i) {
    result += x[i] * x[i];
  }
  for (size_t i = nb_long_axes; i < number_of_variables_; ++i) {
    result += condition * x[i] * x[i];
  }
  return result;
}

/**
 * @brief The bent cigar function f12 of the bbob-largescale suite, generalized to n / 40 long axes and with a
 * permuted block rotation.
 *
 * f(x) = f_bent_cigar_generalized(R T_asy^0.5(R (x - xopt))) / n + fopt with the permuted block rotation R.
 */
class BbobBentCigarGeneralizedPermblockdiag final : public BbobProblem<BbobBentCigarGeneralizedPermblockdiag> {
public:
  BbobBentCigarGeneralizedPermblockdiag(const size_t function, const size_t dimension, const size_t instance,
                                        const long rseed, const BbobVariant variant, const char *problem_id_template,
                                        const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "3-ill-conditioned", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed + 1000000, dimension, variant)),
        rotation_(dimension, rseed + 2000000, rseed + 1000000, rseed + 3000000), asymmetric_(0.5, dimension),
        bent_cigar_(dimension, 40), fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension), z2_(dimension) {
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
    return bent_cigar_.Raw(z1_) * NormByDimFactor(number_of_variables_) + fopt_;
  }

private:
  const VarsShift shift_;
  const VarsPermblockdiag rotation_;
  const VarsAsymmetric asymmetric_;
  const BentCigarGeneralized bent_cigar_;
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

} // namespace coco
