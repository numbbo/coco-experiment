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
 * @brief Data of the large-scale Gallagher function, shared with transform_vars_gallagher_blockrotation.
 */
struct f_gallagher_shared_data_t {
  size_t number_of_peaks = 0;
  std::vector<coco_problem_ptr_t> sub_problems; /**< @brief Owned by the large-scale Gallagher problem */
  std::vector<double> rotated_x;
  std::vector<size_t> block_sizes, block_size_map, first_non_zero_map;
  coco_matrix_t B;
};

/**
 * @brief Computes rotated_x in the shared data of the large-scale Gallagher function.
 */
struct transform_vars_gallagher_blockrotation_t final : public coco_problem_transformed_t {
  std::vector<double> x_copy;
  f_gallagher_shared_data_t *const shared; /**< @brief Owned by the Gallagher function */

  transform_vars_gallagher_blockrotation_t(coco_problem_ptr_t inner, f_gallagher_shared_data_t *shared)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_gallagher_blockrotation"),
        x_copy(inner_problem->number_of_variables), shared(shared) {}

  /**
   * @brief Evaluates the transformation.
   */
  void evaluate_function(const double *x, double *y) override {
    size_t i, j;

    for (i = 0; i < inner_problem->number_of_variables; ++i) {
      shared->rotated_x[i] = 0;
      for (j = shared->first_non_zero_map[i]; j < shared->first_non_zero_map[i] + shared->block_size_map[i]; ++j) {
        shared->rotated_x[i] += shared->B[i][j - shared->first_non_zero_map[i]] * x[j];
      }
      /*shared->rotated_x[i] = x[i];*/
      x_copy[i] = x[i]; /* to avoid pointer problems*/
    }
    coco_evaluate_function(
        inner_problem.get(), x_copy.data(),
        y); /* does not modify the argument of the call since rotated_x will be used later in the sub_problems*/
    /* this function serves only to compute rotated_x on the problem level, not for each sub-problem
     */
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

static COCO_UNUSED coco_problem_ptr_t transform_vars_gallagher_blockrotation(coco_problem_ptr_t inner_problem,
                                                                             f_gallagher_shared_data_t *shared);
static coco_problem_ptr_t transform_vars_gallagher_blockrotation(coco_problem_ptr_t inner_problem,
                                                                 f_gallagher_shared_data_t *shared) {
  if (coco_problem_best_parameter_not_zero(inner_problem.get())) {
    coco_debug("transform_vars_gallagher_blockrotation(): 'best_parameter' not updated, set to NAN");
    coco_vector_set_to_nan(inner_problem->best_parameter.data(), inner_problem->number_of_variables);
  }
  return std::make_unique<transform_vars_gallagher_blockrotation_t>(std::move(inner_problem), shared);
}
