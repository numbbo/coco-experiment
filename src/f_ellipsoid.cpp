/**
 * @file f_ellipsoid.cpp
 * @brief Implementation of the ellipsoid function and problem.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_affine.cpp"

/**
 * @brief The ellipsoid function.
 */
class f_ellipsoid_t final : public coco_problem_s {
public:
  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double raw(std::span<const double> x) const;

  f_ellipsoid_t(const size_t number_of_variables, const double conditioning)
      : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("ellipsoid function", -5.0, 5.0, 0.0);
    set_id("%s_d%02zu", "ellipsoid", number_of_variables);

    scales.resize(number_of_variables);
    for (size_t i = 0; i < number_of_variables; ++i) {
      const double exponent = 1.0 * static_cast<double>(i) / (static_cast<double>(number_of_variables) - 1.0);
      scales[i] = pow(conditioning, exponent);
    }

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

  /**
   * @brief Evaluates the gradient of the ellipsoid function.
   */
  void evaluate_gradient(std::span<const double> x, std::span<double> y) override {
    for (size_t i = 0; i < number_of_variables; ++i)
      y[i] = 2.0 * scales[i] * x[i];
  }

private:
  std::vector<double> scales; /**< @brief conditioning^(i / (n - 1)) */
};

/**
 * @brief Implements the ellipsoid function without connections to any COCO
 * structures.
 */
double f_ellipsoid_t::raw(std::span<const double> x) const {

  double result = x[0] * x[0];
  for (size_t i = 1; i < number_of_variables; ++i)
    result += scales[i] * x[i] * x[i];

  return result;
}

/**
 * @brief The BBOB ellipsoid function f2, also in the bbob-boxed and bbob-largescale suites.
 *
 * f(x) = f_ellipsoid^1e6(T_osz(x - xopt)) + fopt, normalized by the dimension before fopt is added in the
 * large-scale suite.
 */
class bbob_ellipsoid_t final : public bbob_problem_t<bbob_ellipsoid_t> {
public:
  bbob_ellipsoid_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                   const coco_bbob_variant_t variant, const char *problem_id_template,
                   const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "1-separable", problem_id_template, problem_name_template),
        shift(bbob_xopt(rseed, dimension, variant)), ellipsoid(dimension, 1.0e6),
        fopt(bbob2009_compute_fopt(function, instance)), normalize(variant == coco_bbob_variant_t::largescale),
        z(dimension) {
    /* The optimum of the ellipsoid is the origin, map it back through the transformations */
    oscillate.inverse(best_parameter, z);
    shift.inverse(z, best_parameter);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    shift.apply(x, z);
    oscillate.apply(z, z);
    double y = ellipsoid.raw(z);
    if (normalize)
      y *= coco_norm_by_dim_factor(number_of_variables);
    return y + fopt;
  }

private:
  const vars_shift_t shift;
  const vars_oscillate_t oscillate{};
  const f_ellipsoid_t ellipsoid;
  const double fopt;
  const bool normalize;
  std::vector<double> z; /**< @brief The transformed decision values */
};

/**
 * @brief The BBOB rotated ellipsoid function f10, also in the bbob-boxed and bbob-noisy suites.
 *
 * f(x) = f_ellipsoid^conditioning(T_osz(R (x - xopt))) + fopt.
 */
class bbob_ellipsoid_rotated_t final : public bbob_problem_t<bbob_ellipsoid_rotated_t> {
public:
  bbob_ellipsoid_rotated_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                           const coco_bbob_variant_t variant, const double conditioning,
                           const char *problem_id_template, const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "3-ill-conditioned", problem_id_template, problem_name_template),
        shift(bbob_xopt(rseed, dimension, variant)),
        rotation(bbob2009_compute_rotation(rseed + 1000000, dimension), std::vector<double>(dimension)),
        ellipsoid(dimension, conditioning), fopt(bbob2009_compute_fopt(function, instance)), z1(dimension),
        z2(dimension) {
    /* The optimum of the ellipsoid is the origin, map it back through the transformations */
    oscillate.inverse(best_parameter, z1);
    rotation.inverse(z1, z2);
    shift.inverse(z2, best_parameter);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    shift.apply(x, z1);
    rotation.apply(z1, z2);
    oscillate.apply(z2, z2);
    return ellipsoid.raw(z2) + fopt;
  }

private:
  const vars_shift_t shift;
  const vars_affine_t rotation;
  const vars_oscillate_t oscillate{};
  const f_ellipsoid_t ellipsoid;
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};

/**
 * @brief The rotated ellipsoid function f10 of the bbob-largescale suite, with a permuted block rotation.
 *
 * f(x) = f_ellipsoid^1e6(T_osz(R (x - xopt))) / n + fopt with the permuted block rotation R.
 */
class bbob_ellipsoid_permblockdiag_t final : public bbob_problem_t<bbob_ellipsoid_permblockdiag_t> {
public:
  bbob_ellipsoid_permblockdiag_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                                 const coco_bbob_variant_t variant, const char *problem_id_template,
                                 const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "3-ill-conditioned", problem_id_template, problem_name_template),
        shift(bbob_xopt(rseed, dimension, variant)),
        rotation(dimension, rseed + 2000000, rseed + 1000000, rseed + 3000000), ellipsoid(dimension, 1.0e6),
        fopt(bbob2009_compute_fopt(function, instance)), z1(dimension), z2(dimension) {
    /* The optimum of the ellipsoid is the origin, map it back through the transformations */
    oscillate.inverse(best_parameter, z1);
    rotation.inverse(z1, z2);
    shift.inverse(z2, best_parameter);
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    shift.apply(x, z1);
    rotation.apply(z1, z2);
    oscillate.apply(z2, z2);
    return ellipsoid.raw(z2) * coco_norm_by_dim_factor(number_of_variables) + fopt;
  }

private:
  const vars_shift_t shift;
  const vars_permblockdiag_t rotation;
  const vars_oscillate_t oscillate{};
  const f_ellipsoid_t ellipsoid;
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};

/**
 * @brief Creates the ellipsoid problem for the constrained BBOB suite
 */
static coco_problem_ptr_t f_ellipsoid_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                 const size_t instance, const long rseed,
                                                                 const coco_bbob_variant_t variant,
                                                                 const char *problem_id_template,
                                                                 const char *problem_name_template) {
  coco_problem_ptr_t problem;

  std::vector<double> xopt(dimension);
  double fopt = bbob2009_compute_fopt(function, instance);
  bbob_compute_xopt(xopt, rseed, variant);

  problem = std::make_unique<f_ellipsoid_t>(dimension, 1.0e6);

  /* TODO (NH): fopt -= problem->evaluate(all_zeros(dimension)) */
  problem = transform_vars_shift(std::move(problem), xopt, 0);
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("1-separable");

  return problem;
}

/**
 * @brief Creates the rotated ellipsoid problem for the constrained
 *        BBOB suite
 */
static coco_problem_ptr_t f_ellipsoid_rotated_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                         const size_t instance, const long rseed,
                                                                         const coco_bbob_variant_t variant,
                                                                         const char *problem_id_template,
                                                                         const char *problem_name_template) {
  coco_problem_ptr_t problem;

  coco_matrix_t M(dimension, dimension);
  std::vector<double> b(dimension);

  std::vector<double> xopt(dimension);
  bbob_compute_xopt(xopt, rseed, variant);
  double fopt = bbob2009_compute_fopt(function, instance);

  const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
  M = rot1;

  problem = std::make_unique<f_ellipsoid_t>(dimension, 1.0e6);
  problem = transform_vars_affine(std::move(problem), M, b);
  problem = transform_vars_shift(std::move(problem), xopt, 0);
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("3-ill-conditioned");

  return problem;
}
