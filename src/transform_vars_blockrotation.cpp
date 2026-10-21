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
#include "transform_functions.cpp"
#include "transform_vars_blockrotation_helpers.cpp"

/**
 * @brief Applies a block-diagonal rotation matrix to the decision values.
 */
class transform_vars_blockrotation_t final : public coco_problem_transformed_t {
public:
  transform_vars_blockrotation_t(coco_problem_ptr_t inner, coco_block_diagonal_matrix_t B,
                                 const bool keep_best_parameter)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_blockrotation"), rotation(std::move(B)),
        Bx(inner_problem->get_number_of_variables()) {
    if (!keep_best_parameter)
      update_best_parameter();

    if (number_of_variables < 100) {
      /* 1e-11 still passes and 1e-12 fails under macOS */
      test(1e-5);
    }
  }

  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    rotation.apply(x, Bx);
    inner_problem->evaluate_function(Bx, y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /*
   * @brief test blockrotation on its own rows and raise coco_error in case
   */
  void test(double precision) const {
    const coco_block_diagonal_matrix_t &B = rotation.matrix();
    std::vector<double> row(number_of_variables), y(number_of_variables);
    for (size_t i = 0; i < number_of_variables; ++i) { /* for each row */
      for (size_t j = 0; j < number_of_variables; ++j)
        row[j] = B(i, j);
      B.multiply(row, y);
      for (size_t j = 0; j < number_of_variables; ++j) { /* check result */
        if (!coco_double_almost_equal(y[j], (i == j) ? 1.0 : 0.0, precision)) {
          coco_error("transform_vars_blockrotation_t::test() with precision %e "
                     "failed on row %zu",
                     precision, i);
        }
      }
    }
  }

  /** @brief Multiplies by the transpose of the (orthogonal) rotation. */
  bool inverse(std::span<const double> inner_x, std::span<double> x) const override {
    rotation.inverse(inner_x, x);
    return true;
  }

private:
  const vars_blockrotation_t rotation;
  std::vector<double> Bx;
};

/**
 * @brief Creates the transformation.
 *
 * If keep_best_parameter is set, the best parameter of the inner problem is kept instead of being computed with
 * inverse(). Only the large-scale Gallagher function uses this, see f_gallagher_ls_t::set_optimum_at_first_peak().
 */
static coco_problem_ptr_t transform_vars_blockrotation(coco_problem_ptr_t inner_problem, coco_block_diagonal_matrix_t B,
                                                       const bool keep_best_parameter = false) {
  assert(B.dimension() == inner_problem->get_number_of_variables());
  return std::make_unique<transform_vars_blockrotation_t>(std::move(inner_problem), std::move(B), keep_best_parameter);
}
