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
  coco_matrix_t B; /**< @brief the block-diagonal matrices*/
  std::vector<double> Bx;
  const size_t dimension;
  std::vector<size_t> block_size_map;    /**< @brief maps a row to the block-size of the block to
                                           which it belong, keep until better way is found */
  std::vector<size_t> first_non_zero_map; /**< @brief maps a row to the index of its first
                                           non zero element */

  transform_vars_blockrotation_t(coco_problem_ptr_t inner, const double *const *B, const size_t number_of_variables,
                                 const size_t *block_sizes, const size_t nb_blocks)
      : coco_problem_transformed_t(std::move(inner), "transform_vars_blockrotation"),
        B(coco_copy_block_matrix(B, number_of_variables, block_sizes, nb_blocks)),
        Bx(inner_problem->number_of_variables), dimension(number_of_variables), block_size_map(number_of_variables),
        first_non_zero_map(number_of_variables) {
    size_t idx_blocksize, next_bs_change, current_blocksize;
    size_t i;

    idx_blocksize = 0;
    next_bs_change = block_sizes[idx_blocksize];
    for (i = 0; i < number_of_variables; i++) {
      if (i >= next_bs_change) {
        idx_blocksize++;
        next_bs_change += block_sizes[idx_blocksize];
      }
      current_blocksize = block_sizes[idx_blocksize];
      block_size_map[i] = current_blocksize;
      first_non_zero_map[i] = next_bs_change - current_blocksize; /* next_bs_change serves also as a
                                                                     cumsum for blocksizes*/
    }

    if (coco_problem_best_parameter_not_zero(inner_problem.get())) {
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

  /*
   * @brief return i-th row of B in y.
   */
  void get_row(size_t i, double *y) const {
    size_t j, current_blocksize, first_non_zero_ind; /* cp-paste from apply */
    current_blocksize = block_size_map[i];
    first_non_zero_ind = first_non_zero_map[i];
    for (j = 0; j < dimension; ++j) {
      y[j] = (j < first_non_zero_ind || j >= first_non_zero_ind + current_blocksize)
                 ? 0
                 : B[i][j - first_non_zero_ind]; /*all B lines start at 0*/
    }
  }

  /*
   * @brief Computes y = Bx.
   */
  void apply(const double *x, double *y) {
    size_t i, j, current_blocksize, first_non_zero_ind;
    assert(x != Bx.data());
    for (i = 0; i < dimension; ++i) {
      current_blocksize = block_size_map[i];
      first_non_zero_ind = first_non_zero_map[i];
      Bx[i] = 0;
      /*compute y[i] = < B[i,:] , x >  */
      for (j = first_non_zero_ind; j < first_non_zero_ind + current_blocksize; ++j) {
        Bx[i] += B[i][j - first_non_zero_ind] * x[j]; /*all B lines start at 0*/
      }
    }
    if (y != Bx.data()) {
      for (i = 0; i < dimension; ++i) {
        y[i] = Bx[i];
      }
    }
  }

  void evaluate_function(const double *x, double *y) override {
    apply(x, Bx.data());
    coco_evaluate_function(inner_problem.get(), Bx.data(), y);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /*
   * @brief test blockrotation on its own rows and raise coco_error in case
   */
  void test(double precision) {
    size_t i, j;
    size_t number_of_variables = coco_problem_get_dimension(this);
    std::vector<double> y(number_of_variables);
    for (i = 0; i < number_of_variables; ++i) { /* for each row */
      get_row(i, y.data());
      apply(y.data(), y.data());
      for (j = 0; j < number_of_variables; ++j) { /* check result */
        if (!coco_double_almost_equal(y[j], (i == j) ? 1.0 : 0.0, precision)) {
          coco_error("transform_vars_blockrotation_test() with precision %e "
                     "failed on row %lu",
                     precision, (unsigned long)i);
        }
        /* printf("%f ", y[j]); */
      }
      /* printf("\n"); */
    }
  }
};

/*
 * @brief Computes y = Bx for the given block-rotation problem.
 */
static void transform_vars_blockrotation_apply(coco_problem_t *problem, const double *x, double *y) {
  assert(dynamic_cast<transform_vars_blockrotation_t *>(problem) != nullptr);
  static_cast<transform_vars_blockrotation_t *>(problem)->apply(x, y);
}

static coco_problem_ptr_t transform_vars_blockrotation(coco_problem_ptr_t inner_problem, const double *const *B,
                                                       const size_t number_of_variables, const size_t *block_sizes,
                                                       const size_t nb_blocks) {
  assert(number_of_variables > 0); /*tmp*/
  return std::make_unique<transform_vars_blockrotation_t>(std::move(inner_problem), B, number_of_variables, block_sizes,
                                                          nb_blocks);
}
