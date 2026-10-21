/**
 * @file f_sharp_ridge.cpp
 * @brief Implementation of the sharp ridge function and problem.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

/**
 * @brief Implements the sharp ridge function without connections to any COCO
 * structures.
 */
static double f_sharp_ridge_raw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  static const double alpha = 100.0;
  const double d_vars_40 = 1.0; /* generalized: number_of_variables <= 40 ? 1 :
                                   number_of_variables / 40.0; */
  const size_t vars_40 = coco_double_to_size_t(ceil(d_vars_40));
  size_t i = 0;

  assert(number_of_variables > 1);

  double result = 0.0;
  for (i = vars_40; i < number_of_variables; ++i) {
    result += x[i] * x[i];
  }
  result = alpha * sqrt(result / d_vars_40);
  for (i = 0; i < vars_40; ++i)
    result += x[i] * x[i] / d_vars_40;

  return result;
}

/**
 * @brief The BBOB sharp ridge function f13, also in the bbob-boxed suite.
 *
 * f(x) = f_sharp_ridge(Q Λ^10 R (x - xopt)) + fopt.
 */
class bbob_sharp_ridge_t final : public bbob_problem_t<bbob_sharp_ridge_t> {
public:
  bbob_sharp_ridge_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                     const coco_bbob_variant_t variant, const char *problem_id_template,
                     const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "3-ill-conditioned", problem_id_template, problem_name_template),
        shift(bbob_xopt(rseed, dimension, variant)),
        conditioned_rotation(conditioned_rotation_matrix(rseed, dimension), std::vector<double>(dimension)),
        fopt(bbob2009_compute_fopt(function, instance)), z1(dimension), z2(dimension) {
    /* The optimum of the sharp ridge function is the origin, map it back through the transformations */
    conditioned_rotation.inverse(best_parameter, z1);
    shift.inverse(z1, best_parameter);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    shift.apply(x, z1);
    conditioned_rotation.apply(z1, z2);
    return f_sharp_ridge_raw(z2) + fopt;
  }

private:
  /** @brief Returns Q Λ^10 R with the rotations Q and R. */
  static coco_matrix_t conditioned_rotation_matrix(const long rseed, const size_t dimension) {
    const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
    const coco_matrix_t rot2 = bbob2009_compute_rotation(rseed, dimension);
    std::vector<double> scales(dimension);
    for (size_t k = 0; k < dimension; ++k) {
      const double exponent = 1.0 * static_cast<int>(k) / (static_cast<double>(dimension) - 1.0);
      scales[k] = pow(sqrt(10), exponent);
    }
    return coco_matrix_product_with_diagonal(rot1, scales, rot2);
  }

  const vars_shift_t shift;
  const vars_affine_t conditioned_rotation;
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};
