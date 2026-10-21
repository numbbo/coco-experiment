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
struct transform_obj_norm_by_dim_t final : public coco_problem_transformed_t {
  explicit transform_obj_norm_by_dim_t(coco_problem_t *inner_problem)
      : coco_problem_transformed_t(inner_problem, "transform_obj_norm_by_dim") {}

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(const double *x, double *y) override {
    coco_evaluate_function(inner_problem, x, y);
    y[0] *= bbob2009_fmin(1, 40. / ((double)number_of_variables));
    /* Wassim: might want to use a function (with no 40) here that we can put in a helpers file */

    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Creates the transformation.
 */
static coco_problem_t *transform_obj_norm_by_dim(coco_problem_t *inner_problem) {
  return new transform_obj_norm_by_dim_t(inner_problem);
}
