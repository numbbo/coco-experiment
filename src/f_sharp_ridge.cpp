/**
 * @file f_sharp_ridge.cpp
 * @brief Implementation of the sharp ridge function and problem.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

namespace coco {

/**
 * @brief Implements the sharp ridge function without connections to any COCO
 * structures.
 */
static double SharpRidgeRaw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  static const double alpha = 100.0;
  const double d_vars_40 = 1.0; /* generalized: number_of_variables <= 40 ? 1 :
                                   number_of_variables / 40.0; */
  const size_t vars_40 = DoubleToSizeT(ceil(d_vars_40));
  size_t i = 0;

  assert(number_of_variables > 1);

  double result = 0.0;
  for (i = vars_40; i < number_of_variables; ++i) {
    result += x[i] * x[i];
  }
  result = alpha * sqrt(result / d_vars_40);
  for (i = 0; i < vars_40; ++i)
    result += x[i] * x[i] / d_vars_40;

  return result;
}

/**
 * @brief The BBOB sharp ridge function f13, also in the bbob-boxed suite.
 *
 * f(x) = f_sharp_ridge(Q Λ^10 R (x - xopt)) + fopt.
 */
class BbobSharpRidge final : public BbobProblem<BbobSharpRidge> {
public:
  BbobSharpRidge(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                 const BbobVariant variant, const char *problem_id_template, const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "3-ill-conditioned", problem_id_template, problem_name_template),
        shift_(BbobXopt(rseed, dimension, variant)),
        conditioned_rotation_(ConditionedRotationMatrix(rseed, dimension), std::vector<double>(dimension)),
        fopt_(Bbob2009ComputeFopt(function, instance)), z1_(dimension), z2_(dimension) {
    /* The optimum of the sharp ridge function is the origin, map it back through the transformations */
    conditioned_rotation_.Inverse(best_parameter_, z1_);
    shift_.Inverse(z1_, best_parameter_);
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    shift_.Apply(x, z1_);
    conditioned_rotation_.Apply(z1_, z2_);
    return SharpRidgeRaw(z2_) + fopt_;
  }

private:
  /** @brief Returns Q Λ^10 R with the rotations Q and R. */
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
  const VarsAffine conditioned_rotation_;
  const double fopt_;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

} // namespace coco
