/**
 * @file transform_vars_permutation.cpp
 * @brief Implementation of permuting the decision values.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "transform_functions.cpp"
#include "transform_vars_permutation_helpers.cpp"

/**
 * @brief Permutes the decision values (or applies the inverse permutation).
 */
class transform_vars_permutation_t final : public coco_problem_transformed_t {
public:
  transform_vars_permutation_t(coco_problem_ptr_t inner, const size_t *P, const size_t number_of_variables,
                               const bool inverse)
      : coco_problem_transformed_t(std::move(inner),
                                   inverse ? "transform_vars_inverse_permutation" : "transform_vars_permutation"),
        permutation(P, number_of_variables, inverse), permuted_x(number_of_variables) {
    update_best_parameter();
  }

  void evaluate_function(const double *x, double *y) override {
    permutation.apply(x, permuted_x.data());
    inner_problem->evaluate_function(permuted_x.data(), y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /** @brief Applies the inverse permutation. */
  bool inverse(const double *inner_x, double *x) const override {
    permutation.inverse(inner_x, x);
    return true;
  }

private:
  const vars_permutation_t permutation;
  std::vector<double> permuted_x; /**< @brief The permuted decision values */
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_vars_permutation(coco_problem_ptr_t inner_problem, const size_t *P,
                                                     const size_t number_of_variables) {
  assert(number_of_variables == inner_problem->get_number_of_variables()); /*tmp*/
  return std::make_unique<transform_vars_permutation_t>(std::move(inner_problem), P, number_of_variables, false);
}

/**
 * @brief Creates the inverse transformation.
 */
static coco_problem_ptr_t transform_vars_inverse_permutation(coco_problem_ptr_t inner_problem, const size_t *P,
                                                             const size_t number_of_variables) {
  assert(number_of_variables == inner_problem->get_number_of_variables());
  return std::make_unique<transform_vars_permutation_t>(std::move(inner_problem), P, number_of_variables, true);
}
