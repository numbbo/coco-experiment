/**
 * @file f_sphere.cpp
 * @brief Implementation of the sphere function and problem.
 */

#include <assert.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_shift.cpp"

/**
 * @brief Implements the sphere function without connections to any COCO
 * structures.
 */
static double f_sphere_raw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  size_t i = 0;

  double result = 0.0;
  for (i = 0; i < number_of_variables; ++i) {
    result += x[i] * x[i];
  }

  return result;
}

/**
 * @brief The sphere function.
 */
class f_sphere_t final : public coco_problem_s {
public:
  explicit f_sphere_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("sphere function", -5.0, 5.0, 0.0);
    set_id("%s_d%02zu", "sphere", number_of_variables);

    /* Compute best solution */
    evaluate_function(best_parameter, best_value);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives == 1);
    y[0] = f_sphere_raw(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the gradient of the sphere function.
   */
  void evaluate_gradient(std::span<const double> x, std::span<double> y) override {

    for (size_t i = 0; i < number_of_variables; ++i) {
      y[i] = 2.0 * x[i];
    }
  }
};

/**
 * @brief The BBOB sphere function f1, also in the bbob-boxed, bbob-largescale and bbob-noisy suites.
 *
 * f(x) = f_sphere(x - xopt) + fopt, normalized by the dimension before fopt is added in the large-scale suite.
 */
class bbob_sphere_t final : public bbob_problem_t<bbob_sphere_t> {
public:
  bbob_sphere_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                const coco_bbob_variant_t variant, const char *problem_id_template, const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "1-separable", problem_id_template, problem_name_template),
        shift(bbob_xopt(rseed, dimension, variant)), fopt(bbob2009_compute_fopt(function, instance)),
        normalize(variant == coco_bbob_variant_t::largescale), z(dimension) {
    /* The optimum of the sphere function is the origin, map it back through the transformation */
    shift.inverse(best_parameter, best_parameter);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    shift.apply(x, z);
    double y = f_sphere_raw(z);
    if (normalize)
      y *= coco_norm_by_dim_factor(number_of_variables);
    return y + fopt;
  }

private:
  const vars_shift_t shift;
  const double fopt;
  const bool normalize;
  std::vector<double> z; /**< @brief The transformed decision values */
};

/**
 * @brief Creates the sphere problem for the constrained BBOB suite.
 */
static coco_problem_ptr_t f_sphere_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                              const size_t instance, const long rseed,
                                                              const char *problem_id_template,
                                                              const char *problem_name_template) {
  coco_problem_ptr_t problem;

  std::vector<double> xopt(dimension);
  bbob_compute_xopt(xopt, rseed, coco_bbob_variant_t::bbob);
  double fopt = bbob2009_compute_fopt(function, instance);

  problem = std::make_unique<f_sphere_t>(dimension);
  problem = transform_vars_shift(std::move(problem), xopt, 0);
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("1-separable");

  return problem;
}
