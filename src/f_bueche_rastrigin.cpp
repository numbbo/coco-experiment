/**
 * @file f_bueche_rastrigin.cpp
 * @brief Implementation of the Bueche-Rastrigin function and problem.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

/**
 * @brief Implements the Bueche-Rastrigin function without connections to any
 * COCO structures.
 */
static double f_bueche_rastrigin_raw(std::span<const double> x) {
  const size_t number_of_variables = x.size();
  double tmp = 0., tmp2 = 0.;

  for (size_t i = 0; i < number_of_variables; ++i) {
    tmp += cos(2 * coco_pi * x[i]);
    tmp2 += x[i] * x[i];
  }
  double result = 10.0 * (static_cast<double>(number_of_variables) - tmp) + tmp2 + 0;
  return result;
}

/**
 * @brief The Bueche-Rastrigin function.
 */
class f_bueche_rastrigin_t final : public coco_problem_s {
public:
  explicit f_bueche_rastrigin_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("Bueche-Rastrigin function", -5.0, 5.0, 0.0);
    set_id("%s_d%02zu", "bueche-rastrigin", number_of_variables);

    /* Compute best solution */
    evaluate_function(best_parameter, best_value);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives == 1);
    y[0] = f_bueche_rastrigin_raw(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief The BBOB Bueche-Rastrigin function f4, also in the bbob-boxed and bbob-largescale suites.
 *
 * f(x) = f_rastrigin(s T_osz(x - xopt)) + fopt + 100 f_pen(x) with the 's_i scaling' s, normalized by the
 * dimension before fopt is added in the large-scale suite.
 */
class bbob_bueche_rastrigin_t final : public bbob_problem_t<bbob_bueche_rastrigin_t> {
public:
  bbob_bueche_rastrigin_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                          const coco_bbob_variant_t variant, const char *problem_id_template,
                          const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "1-separable", problem_id_template, problem_name_template),
        shift(xopt(rseed, dimension, variant)), brs(dimension), fopt(bbob2009_compute_fopt(function, instance)),
        normalize(variant == coco_bbob_variant_t::largescale), z(dimension) {
    /* The optimum of the Bueche-Rastrigin function is the origin, map it back through the transformations */
    brs.inverse(best_parameter, z);
    oscillate.inverse(z, z);
    shift.inverse(z, best_parameter);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = coco_boundary_penalty(x, smallest_values_of_interest, largest_values_of_interest);
    shift.apply(x, z);
    oscillate.apply(z, z);
    brs.apply(z, z);
    double y = f_bueche_rastrigin_raw(z);
    if (normalize)
      y *= coco_norm_by_dim_factor(number_of_variables);
    y += fopt;
    return y + penalty_factor * penalty;
  }

private:
  /** @brief The location of the optimum: xopt with the absolute value of every other component. */
  static std::vector<double> xopt(const long rseed, const size_t dimension, const coco_bbob_variant_t variant) {
    std::vector<double> xopt = bbob_xopt(rseed, dimension, variant);
    /* OME: This step is in the legacy C code but _not_ in the function
     * description. */
    for (size_t i = 0; i < dimension; i += 2) {
      xopt[i] = fabs(xopt[i]);
    }
    return xopt;
  }

  static constexpr double penalty_factor = 100.0;

  const vars_shift_t shift;
  const vars_oscillate_t oscillate{};
  const vars_brs_t brs;
  const double fopt;
  const bool normalize;
  std::vector<double> z; /**< @brief The transformed decision values */
};
