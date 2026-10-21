/**
 * @file f_bent_cigar_generalized.cpp
 * @brief Implementation of the generalized bent cigar function and problem.
 */

#include <assert.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

/**
 * @brief The generalized bent cigar function.
 */
class f_bent_cigar_generalized_t final : public coco_problem_s {
public:
  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double raw(std::span<const double> x) const;

  f_bent_cigar_generalized_t(const size_t number_of_variables, size_t proportion_long_axes_denom)
      : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("generalized bent cigar function", -5.0, 5.0, 0.0);
    set_id("%s_d%04zu", "bent_cigar", number_of_variables);
    this->proportion_long_axes_denom = proportion_long_axes_denom;

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
  size_t proportion_long_axes_denom{};
};

/**
 * @brief Implements the generalized bent cigar function without connections to
 * any COCO structures.
 */
double f_bent_cigar_generalized_t::raw(std::span<const double> x) const {

  static const double condition = 1.0e6;
  double result = 0;
  size_t nb_long_axes = number_of_variables / proportion_long_axes_denom;
  if (number_of_variables % proportion_long_axes_denom != 0) {
    nb_long_axes += 1;
  }
  for (size_t i = 0; i < nb_long_axes; ++i) {
    result += x[i] * x[i];
  }
  for (size_t i = nb_long_axes; i < number_of_variables; ++i) {
    result += condition * x[i] * x[i];
  }
  return result;
}

/**
 * @brief The bent cigar function f12 of the bbob-largescale suite, generalized to n / 40 long axes and with a
 * permuted block rotation.
 *
 * f(x) = f_bent_cigar_generalized(R T_asy^0.5(R (x - xopt))) / n + fopt with the permuted block rotation R.
 */
class bbob_bent_cigar_generalized_permblockdiag_t final
    : public bbob_problem_t<bbob_bent_cigar_generalized_permblockdiag_t> {
public:
  bbob_bent_cigar_generalized_permblockdiag_t(const size_t function, const size_t dimension, const size_t instance,
                                              const long rseed, const coco_bbob_variant_t variant,
                                              const char *problem_id_template, const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "3-ill-conditioned", problem_id_template, problem_name_template),
        shift(bbob_xopt(rseed + 1000000, dimension, variant)),
        rotation(dimension, rseed + 2000000, rseed + 1000000, rseed + 3000000), asymmetric(0.5, dimension),
        bent_cigar(dimension, 40), fopt(bbob2009_compute_fopt(function, instance)), z1(dimension), z2(dimension) {
    /* The optimum of the bent cigar function is the origin, map it back through the transformations */
    rotation.inverse(best_parameter, z1);
    asymmetric.inverse(z1, z1);
    rotation.inverse(z1, z2);
    shift.inverse(z2, best_parameter);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    shift.apply(x, z1);
    rotation.apply(z1, z2);
    asymmetric.apply(z2, z2);
    rotation.apply(z2, z1);
    return bent_cigar.raw(z1) * coco_norm_by_dim_factor(number_of_variables) + fopt;
  }

private:
  const vars_shift_t shift;
  const vars_permblockdiag_t rotation;
  const vars_asymmetric_t asymmetric;
  const f_bent_cigar_generalized_t bent_cigar;
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};
