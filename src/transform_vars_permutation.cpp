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
  transform_vars_permutation_t(coco_problem_ptr_t inner, std::span<const size_t> P, const bool inverse)
      : coco_problem_transformed_t(std::move(inner),
                                   inverse ? "transform_vars_inverse_permutation" : "transform_vars_permutation"),
        permutation(P, inverse), permuted_x(P.size()) {
    update_best_parameter();
  }

  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    permutation.apply(x, permuted_x);
    inner_problem->evaluate_function(permuted_x, y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /** @brief Applies the inverse permutation. */
  bool inverse(std::span<const double> inner_x, std::span<double> x) const override {
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
static coco_problem_ptr_t transform_vars_permutation(coco_problem_ptr_t inner_problem, std::span<const size_t> P) {
  assert(P.size() == inner_problem->get_number_of_variables());
  return std::make_unique<transform_vars_permutation_t>(std::move(inner_problem), P, false);
}

/**
 * @brief Creates the inverse transformation.
 */
static coco_problem_ptr_t transform_vars_inverse_permutation(coco_problem_ptr_t inner_problem,
                                                             std::span<const size_t> P) {
  assert(P.size() == inner_problem->get_number_of_variables());
  return std::make_unique<transform_vars_permutation_t>(std::move(inner_problem), P, true);
}
