/**
 * @file transform_vars_affine.cpp
 * @brief Implementation of performing an affine transformation on decision values.
 *
 * x |-> Mx + b
 *
 * Currently, the best parameter is transformed correctly only in the simple
 * cases where M is orthogonal which is always the case for the `bbob`
 * functions. How to code this for general transformations of the above form,
 * see https://github.com/numbbo/coco/issues/814#issuecomment-303724400
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "transform_functions.cpp"

namespace coco {

/**
 * @brief Applies the affine transformation x |-> Mx + b to the decision values.
 */
class TransformVarsAffine final : public TransformedProblem {
public:
  TransformVarsAffine(ProblemPtr inner, Matrix M, std::vector<double> b)
      : TransformedProblem(std::move(inner), "transform_vars_affine"), affine_(std::move(M), std::move(b)),
        transformed_x_(inner_problem_->number_of_variables()) {
    /*
     * TODOs:
     * - Calculate new smallest/largest values of interest?
     * - Resize bounds vectors if input and output dimensions do not match
     */
    assert(affine_.matrix().number_of_rows() == inner_problem_->number_of_variables());
    UpdateBestParameter();
  }

  /**
   * @brief Evaluates the transformed objective function.
   */
  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    if (inner_problem_->number_of_objectives() == 0) {
      inner_problem_->EvaluateFunction(x, y);
      return;
    }

    affine_.Apply(x, transformed_x_);
    inner_problem_->EvaluateFunction(transformed_x_, y);

    AssertNotBelowBestValue(x, y);
  }

  /**
   * @brief Evaluates the transformed constraint.
   */
  void EvaluateConstraint(std::span<const double> x, std::span<double> y, int update_counter) override {
    if (inner_problem_->number_of_constraints() == 0) {
      inner_problem_->EvaluateConstraint(x, y, update_counter);
      return;
    }

    affine_.Apply(x, transformed_x_);
    inner_problem_->EvaluateConstraint(transformed_x_, y, update_counter);
  }

  /**
   * @brief Evaluates the gradient of the transformed function.
   */
  void EvaluateGradient(std::span<const double> x, std::span<double> y) override {

    std::vector<double> gradient(inner_problem_->number_of_variables(), 0.0);

    affine_.Apply(x, transformed_x_);
    inner_problem_->EvaluateGradient(transformed_x_, y);

    /* grad_(f o g )(x), where g(x) = M * x + b, equals to
     * M^T * grad_f(M *x + b)
     */
    const Matrix &M = affine_.matrix();
    for (size_t j = 0; j < inner_problem_->number_of_variables(); ++j) {
      for (size_t i = 0; i < inner_problem_->number_of_variables(); ++i) {
        gradient[j] += y[i] * M(i, j);
      }
    }
    for (size_t i = 0; i < inner_problem_->number_of_variables(); ++i)
      y[i] = gradient[i];
  }

  /**
   * @brief Solves M x = inner_x - b, which is M^T (inner_x - b) if M is orthogonal. Not possible if M is singular.
   */
  bool Inverse(std::span<const double> inner_x, std::span<double> x) const override {
    return affine_.Inverse(inner_x, x);
  }

private:
  const VarsAffine affine_;
  std::vector<double> transformed_x_; /**< @brief The transformed decision values */
};

} // namespace coco
