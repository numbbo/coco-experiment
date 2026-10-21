/**
 * @file transform_obj_norm_by_dim.cpp
 * @brief Implementation of normalizing the raw fitness functions by the dimensions
 * Mostly used to in the large-scale testsuite
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"

/**
 * @brief Normalizes the objective value of the inner problem by the dimension.
 */
class transform_obj_norm_by_dim_t final : public coco_problem_transformed_t {
public:
  explicit transform_obj_norm_by_dim_t(coco_problem_ptr_t inner)
      : coco_problem_transformed_t(std::move(inner), "transform_obj_norm_by_dim") {}

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(const double *x, double *y) override {
    inner_problem->evaluate_function(x, y);
    y[0] *= bbob2009_fmin(1, 40. / (static_cast<double>(number_of_variables)));
    /* Wassim: might want to use a function (with no 40) here that we can put in a helpers file */

    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /** @brief The decision values are not transformed. */
  bool inverse(const double *inner_x, double *x) const override { return identity_inverse(inner_x, x); }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_ptr_t transform_obj_norm_by_dim(coco_problem_ptr_t inner_problem) {
  return std::make_unique<transform_obj_norm_by_dim_t>(std::move(inner_problem));
}
