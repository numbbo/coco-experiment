/**
 * @file f_sharp_ridge_generalized.cpp
 * @brief Implementation of the generalized sharp ridge function and problem.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

/**
 * @brief The sharp ridge function.
 */
class f_sharp_ridge_generalized_t final : public coco_problem_s {
public:
  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double raw(std::span<const double> x) const;

  f_sharp_ridge_generalized_t(const size_t number_of_variables, size_t proportion_of_linear_dims)
      : coco_problem_s(number_of_variables, 1, 0) {
    /* Wassim: proportion_of_linear_dims should probably be allowed to be non-integer */
    init_from_scalars("sharp ridge function", -5.0, 5.0, 0.0);
    set_id("%s_d%02zu", "sharp_ridge_generalized", number_of_variables);
    this->proportion_of_linear_dims = proportion_of_linear_dims;

    /* Compute best solution */
    evaluate_function(best_parameter, best_value);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives == 1);
    y[0] = raw(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

private:
  size_t proportion_of_linear_dims{};
};

/**
 * @brief Implements the generalized sharp ridge function without connections to
 * any COCO structures.
 */
double f_sharp_ridge_generalized_t::raw(std::span<const double> x) const {

  static const double alpha = 100.0;
  size_t i = 0, number_linear_dimensions;

  assert(number_of_variables > 1);

  double result = 0.0;
  number_linear_dimensions = number_of_variables / proportion_of_linear_dims;
  if (number_of_variables % proportion_of_linear_dims != 0) {
    number_linear_dimensions += 1;
  }
  for (i = number_linear_dimensions; i < number_of_variables; ++i) {
    result += x[i] * x[i];
  }
  result = alpha * sqrt(result);
  for (i = 0; i < number_linear_dimensions; i++) {
    result += x[i] * x[i];
  }

  return result;
}

/**
 * @brief The sharp ridge function f13 of the bbob-largescale suite, generalized to n / 40 linear dimensions and
 * with permuted block rotations.
 *
 * f(x) = f_sharp_ridge_generalized(Q Λ^10 R (x - xopt)) / n + fopt with the permuted block rotations R and Q.
 */
class bbob_sharp_ridge_generalized_permblockdiag_t final
    : public bbob_problem_t<bbob_sharp_ridge_generalized_permblockdiag_t> {
public:
  bbob_sharp_ridge_generalized_permblockdiag_t(const size_t function, const size_t dimension, const size_t instance,
                                               const long rseed, const coco_bbob_variant_t variant,
                                               const char *problem_id_template, const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "3-ill-conditioned", problem_id_template, problem_name_template),
        shift(bbob_xopt(rseed, dimension, variant)), rotation1(dimension, rseed + 5000000, rseed, rseed + 6000000),
        conditioning(10.0, dimension), rotation2(dimension, rseed + 3000000, rseed + 1000000, rseed + 4000000),
        sharp_ridge(dimension, 40), fopt(bbob2009_compute_fopt(function, instance)), z1(dimension), z2(dimension) {
    /* The optimum of the sharp ridge function is the origin, map it back through the transformations */
    rotation2.inverse(best_parameter, z1);
    conditioning.inverse(z1, z1);
    rotation1.inverse(z1, z2);
    shift.inverse(z2, best_parameter);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    shift.apply(x, z1);
    rotation1.apply(z1, z2);
    conditioning.apply(z2, z2);
    rotation2.apply(z2, z1);
    return sharp_ridge.raw(z1) * coco_norm_by_dim_factor(number_of_variables) + fopt;
  }

private:
  const vars_shift_t shift;
  const vars_permblockdiag_t rotation1;
  const vars_conditioning_t conditioning;
  const vars_permblockdiag_t rotation2;
  const f_sharp_ridge_generalized_t sharp_ridge;
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};
