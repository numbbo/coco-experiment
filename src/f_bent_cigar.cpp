/**
 * @file f_bent_cigar.cpp
 * @brief Implementation of the bent cigar function and problem.
 */

#include <assert.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_affine.cpp"
#include "transform_vars_shift.cpp"

/**
 * @brief Implements the bent cigar function without connections to any COCO
 * structures.
 */
static double f_bent_cigar_raw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  static const double condition = 1.0e6;

  double result = x[0] * x[0];
  for (size_t i = 1; i < number_of_variables; ++i) {
    result += condition * x[i] * x[i];
  }
  return result;
}

/**
 * @brief The bent cigar function.
 */
class f_bent_cigar_t final : public coco_problem_s {
public:
  explicit f_bent_cigar_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("bent cigar function", -5.0, 5.0, 0.0);
    set_id("%s_d%02zu", "bent_cigar", number_of_variables);

    /* Compute best solution */
    evaluate_function(best_parameter, best_value);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives == 1);
    y[0] = f_bent_cigar_raw(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the gradient of the bent cigar function.
   */
  void evaluate_gradient(std::span<const double> x, std::span<double> y) override {

    static const double condition = 1.0e6;

    y[0] = 2.0 * x[0];
    for (size_t i = 1; i < number_of_variables; ++i)
      y[i] = 2.0 * condition * x[i];
  }
};

/**
 * @brief The BBOB bent cigar function f12, also in the bbob-boxed suite.
 *
 * f(x) = f_bent_cigar(R T_asy^0.5(R (x - xopt))) + fopt.
 */
class bbob_bent_cigar_t final : public bbob_problem_t<bbob_bent_cigar_t> {
public:
  bbob_bent_cigar_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                    const coco_bbob_variant_t variant, const char *problem_id_template,
                    const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "3-ill-conditioned", problem_id_template, problem_name_template),
        shift(bbob_xopt(rseed + 1000000, dimension, variant)),
        rotation(bbob2009_compute_rotation(rseed + 1000000, dimension), std::vector<double>(dimension)),
        asymmetric(0.5, dimension), fopt(bbob2009_compute_fopt(function, instance)), z1(dimension), z2(dimension) {
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
    return f_bent_cigar_raw(z1) + fopt;
  }

private:
  const vars_shift_t shift;
  const vars_affine_t rotation;
  const vars_asymmetric_t asymmetric;
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};

/**
 * @brief Creates the bent cigar problem for the constrained BBOB suite.
 */
static coco_problem_ptr_t f_bent_cigar_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                  const size_t instance, const long rseed,
                                                                  const coco_bbob_variant_t variant,
                                                                  const char *problem_id_template,
                                                                  const char *problem_name_template) {

  coco_problem_ptr_t problem;

  coco_matrix_t M(dimension, dimension);
  std::vector<double> b(dimension);

  std::vector<double> xopt(dimension);
  double fopt = bbob2009_compute_fopt(function, instance);
  bbob_compute_xopt(xopt, rseed + 1000000, variant);

  const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
  M = rot1;

  problem = std::make_unique<f_bent_cigar_t>(dimension);
  problem = transform_obj_shift(std::move(problem), fopt);
  problem = transform_vars_affine(std::move(problem), M, b);
  problem = transform_vars_shift(std::move(problem), xopt, 0);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("3-ill-conditioned");

  return problem;
}
