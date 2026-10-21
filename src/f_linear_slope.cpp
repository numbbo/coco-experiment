/**
 * @file f_linear_slope.cpp
 * @brief Implementation of the linear slope function and problem.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_shift.cpp"

/**
 * @brief Implements the linear slope function without connections to any COCO
 * structures.
 */
static double f_linear_slope_raw(std::span<const double> x, std::span<const double> best_parameter,
                                 std::span<const double> slopes) {
  double result = 0.0;

  assert(best_parameter.size() == x.size() && slopes.size() == x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    const double si = best_parameter[i] > 0.0 ? slopes[i] : -slopes[i];
    /* boundary handling */
    if (x[i] * best_parameter[i] < 25.0) {
      result += 5.0 * fabs(si) - si * x[i];
    } else {
      result += 5.0 * fabs(si) - si * best_parameter[i];
    }
  }

  return result;
}

/**
 * @brief The linear slope function.
 */
class f_linear_slope_t final : public coco_problem_s {
public:
  f_linear_slope_t(const size_t number_of_variables, std::span<const double> best_parameter)
      : coco_problem_s(number_of_variables, 1, 0), slopes(number_of_variables) {
    assert(best_parameter.size() == number_of_variables);
    static const double alpha = 100.0;
    for (size_t i = 0; i < number_of_variables; ++i)
      slopes[i] = pow(sqrt(alpha), static_cast<double>(i) / (static_cast<double>(number_of_variables) - 1));

    /* best_parameter will be overwritten below */
    init_from_scalars("linear slope function", -5.0, 5.0, 0.0);
    set_id("%s_d%02zu", "linear_slope", number_of_variables);

    /* Compute best solution */
    for (size_t i = 0; i < number_of_variables; ++i) {
      if (best_parameter[i] < 0.0) {
        this->best_parameter[i] = smallest_values_of_interest[i];
      } else {
        this->best_parameter[i] = largest_values_of_interest[i];
      }
    }
    evaluate_function(this->best_parameter, best_value);
  }

  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double raw(std::span<const double> x) const { return f_linear_slope_raw(x, best_parameter, slopes); }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives == 1);
    y[0] = raw(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the gradient of the linear slope function.
   */
  void evaluate_gradient(std::span<const double>, std::span<double> y) override {
    for (size_t i = 0; i < number_of_variables; ++i)
      y[i] = best_parameter[i] > 0.0 ? -slopes[i] : slopes[i];
  }

private:
  std::vector<double> slopes; /**< @brief |s_i| = sqrt(alpha)^(i / (n - 1)) */
};

/**
 * @brief The BBOB linear slope function f5, also in the bbob-boxed and bbob-largescale suites.
 *
 * f(x) = f_linear_slope(x) + fopt with the optimum at the corner of the region of interest given by the signs of
 * xopt, normalized by the dimension before fopt is added in the large-scale suite.
 */
class bbob_linear_slope_t final : public bbob_problem_t<bbob_linear_slope_t> {
public:
  bbob_linear_slope_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                      const coco_bbob_variant_t variant, const char *problem_id_template,
                      const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "1-separable", problem_id_template, problem_name_template),
        linear_slope(dimension, bbob_xopt(rseed, dimension, variant)), fopt(bbob2009_compute_fopt(function, instance)),
        normalize(variant == coco_bbob_variant_t::largescale) {
    best_parameter = linear_slope.get_best_parameter();
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    double y = linear_slope.raw(x);
    if (normalize)
      y *= coco_norm_by_dim_factor(number_of_variables);
    return y + fopt;
  }

private:
  const f_linear_slope_t linear_slope;
  const double fopt;
  const bool normalize;
};

/**
 * @brief Creates the linear slope problem for the constrained BBOB suite.
 */
static coco_problem_ptr_t f_linear_slope_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                    const size_t instance, const long rseed,
                                                                    const char *problem_id_template,
                                                                    const char *problem_name_template) {
  coco_problem_ptr_t problem;

  std::vector<double> xopt(dimension);
  bbob_compute_xopt(xopt, rseed, coco_bbob_variant_t::bbob);
  double fopt = bbob2009_compute_fopt(function, instance);

  problem = std::make_unique<f_linear_slope_t>(dimension, xopt);
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("1-separable");

  return problem;
}
