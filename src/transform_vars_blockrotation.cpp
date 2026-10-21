/**
 * @file transform_vars_blockrotation.cpp
 * @brief Implementation of performing a block-rotation transformation on
 * decision values.
 *
 * x |-> Bx
 * The matrix B is stored in a 2-D array. Only the content of the blocks are
 * stored
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "transform_vars_blockrotation_helpers.cpp"

/**
 * @brief Applies a block-diagonal rotation matrix to the decision values.
 */
struct transform_vars_blockrotation_t final : public coco_problem_transformed_t {
  const coco_block_diagonal_matrix_t B; /**< @brief The block-diagonal rotation matrix */
  std::vector<double> Bx;

  transform_vars_blockrotation_t(coco_problem_ptr_t inner, coco_block_diagonal_matrix_t B)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_blockrotation"), B(std::move(B)),
        Bx(inner_problem->number_of_variables) {

    if (inner_problem->best_parameter_not_zero()) {
      coco_debug("transform_vars_blockrotation(): 'best_parameter' not updated, "
                 "set to NAN");
      coco_vector_set_to_nan(inner_problem->best_parameter.data(), inner_problem->number_of_variables);
      /* The best parameter is copied from the inner problem */
      best_parameter = inner_problem->best_parameter;
    }

    if (number_of_variables < 100) {
      /* 1e-11 still passes and 1e-12 fails under macOS */
      test(1e-5);
    }
  }

  void evaluate_function(const double *x, double *y) override {
    B.multiply(x, Bx.data());
    coco_evaluate_function(inner_problem.get(), Bx.data(), y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /*
   * @brief test blockrotation on its own rows and raise coco_error in case
   */
  void test(double precision) const {
    std::vector<double> row(number_of_variables), y(number_of_variables);
    for (size_t i = 0; i < number_of_variables; ++i) { /* for each row */
      for (size_t j = 0; j < number_of_variables; ++j)
        row[j] = B(i, j);
      B.multiply(row.data(), y.data());
      for (size_t j = 0; j < number_of_variables; ++j) { /* check result */
        if (!coco_double_almost_equal(y[j], (i == j) ? 1.0 : 0.0, precision)) {
          coco_error("transform_vars_blockrotation_t::test() with precision %e "
                     "failed on row %lu",
                     precision, (unsigned long)i);
        }
      }
    }
  }
};

static coco_problem_ptr_t transform_vars_blockrotation(coco_problem_ptr_t inner_problem,
                                                       coco_block_diagonal_matrix_t B) {
  assert(B.dimension() == inner_problem->number_of_variables);
  return std::make_unique<transform_vars_blockrotation_t>(std::move(inner_problem), std::move(B));
}
