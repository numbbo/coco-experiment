/**
 * @file f_different_powers.cpp
 * @brief Implementation of the different powers function and problem.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_affine.cpp"
#include "transform_vars_shift.cpp"

/**
 * @brief Implements the different powers function without connections to any
 * COCO structures.
 */
static double f_different_powers_raw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  double sum = 0.0;

  for (size_t i = 0; i < number_of_variables; ++i) {
    double exponent = 2.0 + (4.0 * static_cast<double>(i)) / (static_cast<double>(number_of_variables) - 1.0);
    sum += pow(fabs(x[i]), exponent);
  }
  double result = sqrt(sum);

  return result;
}

/**
 * @brief Implements the sign function.
 */
double sign(double x) {

  if (x > 0)
    return 1;
  if (x < 0)
    return -1;
  return 0;
}

/**
 * @brief The different powers function.
 */
class f_different_powers_t final : public coco_problem_s {
public:
  explicit f_different_powers_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("different powers function", -5.0, 5.0, 0.0);
    set_id("%s_d%02zu", "different_powers", number_of_variables);

    /* Compute best solution */
    evaluate_function(best_parameter, best_value);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives == 1);
    y[0] = f_different_powers_raw(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the gradient of the function "different powers".
   */
  void evaluate_gradient(std::span<const double> x, std::span<double> y) override {

    double sum = 0.0;
    double aux;

    for (size_t i = 0; i < number_of_variables; ++i) {
      aux = 2.0 + (4.0 * static_cast<double>(i)) / (static_cast<double>(number_of_variables) - 1.0);
      sum += pow(fabs(x[i]), aux);
    }

    for (size_t i = 0; i < number_of_variables; ++i) {
      aux = 2.0 + (4.0 * static_cast<double>(i)) / (static_cast<double>(number_of_variables) - 1.0);
      y[i] = 0.5 * (aux) / (sum);
      aux -= 1.0;
      y[i] *= pow(fabs(x[i]), aux) * sign(x[i]);
    }
  }
};

/**
 * @brief The BBOB different powers function f14, also in the bbob-boxed and bbob-noisy suites.
 *
 * f(x) = f_different_powers(R (x - xopt)) + fopt.
 */
class bbob_different_powers_t final : public bbob_problem_t<bbob_different_powers_t> {
public:
  bbob_different_powers_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                          const coco_bbob_variant_t variant, const char *problem_id_template,
                          const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "3-ill-conditioned", problem_id_template, problem_name_template),
        shift(bbob_xopt(rseed, dimension, variant)),
        rotation(bbob2009_compute_rotation(rseed + 1000000, dimension), std::vector<double>(dimension)),
        fopt(bbob2009_compute_fopt(function, instance)), z1(dimension), z2(dimension) {
    /* The optimum of the different powers function is the origin, map it back through the transformations */
    rotation.inverse(best_parameter, z1);
    shift.inverse(z1, best_parameter);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    shift.apply(x, z1);
    rotation.apply(z1, z2);
    return f_different_powers_raw(z2) + fopt;
  }

private:
  const vars_shift_t shift;
  const vars_affine_t rotation;
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};

/**
 * @brief Creates the bbob-constrained different powers problem.
 */
static coco_problem_ptr_t f_different_powers_bbob_constrained_problem_allocate(
    const size_t function, const size_t dimension, const size_t instance, const long rseed,
    const coco_bbob_variant_t variant, const char *problem_id_template, const char *problem_name_template) {
  /* Different powers function used in bbob-constrained test suite.
   * In this version, the (unconstrained) optimum, xopt, is set to
   * a distance of 1e-2 to the origin. By doing so, the optimum of
   * the constrained problem is at a "reasonable" distance from
   * the unconstrained one and, hence, the constrained problem is not too easy.
   */

  double fopt, result;
  coco_problem_ptr_t problem;

  coco_matrix_t M(dimension, dimension);
  std::vector<double> b(dimension);

  std::vector<double> xopt(dimension);
  fopt = bbob2009_compute_fopt(function, instance);
  bbob_compute_xopt(xopt, rseed, variant);

  /* Compute Euclidean norm of xopt */
  /* CAVEAT: this implementation is not ideal for large dimensions */
  result = 0.0;
  for (size_t i = 0; i < dimension; ++i) {
    result += xopt[i] * xopt[i];
  }
  result = sqrt(result);

  /* Scale xopt such that the distance to the origin is 1e-2 */
  for (size_t i = 0; i < dimension; ++i) {
    xopt[i] *= 1e-2 / result;
  }

  const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
  M = rot1;

  problem = std::make_unique<f_different_powers_t>(dimension);
  problem = transform_obj_shift(std::move(problem), fopt);
  problem = transform_vars_affine(std::move(problem), M, b);
  problem = transform_vars_shift(std::move(problem), xopt, 0);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("3-ill-conditioned");

  return problem;
}

/**
 * @brief The different powers function f14 of the bbob-largescale suite, with a permuted block rotation.
 *
 * f(x) = f_different_powers(R (x - xopt)) / n + fopt with the permuted block rotation R.
 */
class bbob_different_powers_permblockdiag_t final : public bbob_problem_t<bbob_different_powers_permblockdiag_t> {
public:
  bbob_different_powers_permblockdiag_t(const size_t function, const size_t dimension, const size_t instance,
                                        const long rseed, const coco_bbob_variant_t variant,
                                        const char *problem_id_template, const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "3-ill-conditioned", problem_id_template, problem_name_template),
        shift(bbob_xopt(rseed, dimension, variant)),
        rotation(dimension, rseed + 2000000, rseed + 1000000, rseed + 3000000),
        fopt(bbob2009_compute_fopt(function, instance)), z1(dimension), z2(dimension) {
    /* The optimum of the different powers function is the origin, map it back through the transformations */
    rotation.inverse(best_parameter, z1);
    shift.inverse(z1, best_parameter);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    shift.apply(x, z1);
    rotation.apply(z1, z2);
    return f_different_powers_raw(z2) * coco_norm_by_dim_factor(number_of_variables) + fopt;
  }

private:
  const vars_shift_t shift;
  const vars_permblockdiag_t rotation;
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};
