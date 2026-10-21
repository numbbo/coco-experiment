/**
 * @file transform_vars_gallagher_blockrotation.cpp
 * @brief Implementation of performing a block-rotation transformation on decision values for the Gallagher function.
 * The block-rotation is applied only once per call and the result is stored such that the sub-problems can access it
 *
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "transform_vars_blockrotation_helpers.cpp"
/* #include "coco_utilities.cpp" */

/**
 * @brief Data type in problem->versatile_data of f_gallagher.cpp
 */
typedef struct f_gallagher_versatile_data_s {
  size_t number_of_peaks = 0;
  std::vector<coco_problem_t *> sub_problems; /**< @brief Owned by the large-scale Gallagher problem */
  std::vector<double> rotated_x;
  std::vector<size_t> block_sizes, block_size_map, first_non_zero_map;
  double **B = NULL;
} f_gallagher_versatile_data_t;

/**
 * @brief Computes rotated_x in the versatile data of the large-scale Gallagher function.
 */
struct transform_vars_gallagher_blockrotation_t final : public coco_problem_transformed_t {
  std::vector<double> x_copy;

  explicit transform_vars_gallagher_blockrotation_t(coco_problem_t *inner_problem)
      : coco_problem_transformed_t(inner_problem, "transform_vars_gallagher_blockrotation"),
        x_copy(inner_problem->number_of_variables) {}

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(const double *x, double *y) override {
    size_t i, j;
    f_gallagher_versatile_data_t *versatile = (f_gallagher_versatile_data_t *)versatile_data;

    for (i = 0; i < inner_problem->number_of_variables; ++i) {
      versatile->rotated_x[i] = 0;
      for (j = versatile->first_non_zero_map[i]; j < versatile->first_non_zero_map[i] + versatile->block_size_map[i];
           ++j) {
        versatile->rotated_x[i] += versatile->B[i][j - versatile->first_non_zero_map[i]] * x[j];
      }
      /*((f_gallagher_versatile_data_t *) problem->versatile_data)->rotated_x[i] = x[i];*/
      x_copy[i] = x[i]; /* to avoid pointer problems*/
    }
    coco_evaluate_function(
        inner_problem, x_copy.data(),
        y); /* does not modify the argument of the call since rotated_x will be used later in the sub_problems*/
    /* this function serves only to compute rotated_x on the problem level, not for each sub-problem
     */
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

static COCO_UNUSED coco_problem_t *transform_vars_gallagher_blockrotation(coco_problem_t *inner_problem);
static coco_problem_t *transform_vars_gallagher_blockrotation(coco_problem_t *inner_problem) {
  if (coco_problem_best_parameter_not_zero(inner_problem)) {
    coco_debug("transform_vars_gallagher_blockrotation(): 'best_parameter' not updated, set to NAN");
    coco_vector_set_to_nan(inner_problem->best_parameter.data(), inner_problem->number_of_variables);
  }
  return new transform_vars_gallagher_blockrotation_t(inner_problem);
}
