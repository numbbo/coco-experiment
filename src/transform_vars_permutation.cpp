/**
 * @file transform_vars_permutation.cpp
 * @brief Implementation of permuting the decision values.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "transform_vars_permutation_helpers.cpp"

/**
 * @brief Permutes the decision values (or applies the inverse permutation).
 */
struct transform_vars_permutation_t final : public coco_problem_transformed_t {
  const bool is_inverse;
  std::vector<double> permuted_x;    /**< @brief The permuted decision values */
  std::vector<size_t> P; /**< @brief the permutation matrices*/

  transform_vars_permutation_t(coco_problem_ptr_t inner, const size_t *P, const size_t number_of_variables,
                               const bool inverse)
      : coco_problem_transformed_t(std::move(inner),
                                   inverse ? "transform_vars_inverse_permutation" : "transform_vars_permutation"),
        is_inverse(inverse), permuted_x(number_of_variables), P(P, P + number_of_variables) {}

  void evaluate_function(const double *x, double *y) override {

    if (is_inverse) {
      for (size_t i = 0; i < inner_problem->number_of_variables; ++i) {
        permuted_x[P[i]] = x[i];
      }
    } else {
      for (size_t i = 0; i < inner_problem->number_of_variables; ++i) {
        permuted_x[i] = x[P[i]];
      }
    }
    coco_evaluate_function(inner_problem.get(), permuted_x.data(), y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /** @brief Applies the inverse permutation. */
  bool inverse(const double *inner_x, double *x) const override {
    for (size_t i = 0; i < number_of_variables; ++i) {
      if (is_inverse)
        x[i] = inner_x[P[i]];
      else
        x[P[i]] = inner_x[i];
    }
    return true;
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_permutation(coco_problem_ptr_t inner_problem, const size_t *P,
                                                     const size_t number_of_variables) {
  assert(number_of_variables == inner_problem->number_of_variables); /*tmp*/
  return std::make_unique<transform_vars_permutation_t>(std::move(inner_problem), P, number_of_variables, false);
}

/**
 * @brief Creates the inverse transformation.
 */
static coco_problem_ptr_t transform_vars_inverse_permutation(coco_problem_ptr_t inner_problem, const size_t *P,
                                                             const size_t number_of_variables) {
  assert(number_of_variables == inner_problem->number_of_variables);
  return std::make_unique<transform_vars_permutation_t>(std::move(inner_problem), P, number_of_variables, true);
}
