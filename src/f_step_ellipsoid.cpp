/**
 * @file f_step_ellipsoid.cpp
 * @brief Implementation of the step ellipsoid function and problem.
 *
 * The BBOB step ellipsoid function intertwines the variable and objective
 * transformations in such a way that it is hard to devise a composition of
 * generic transformations to implement it. In the end one would have to
 * implement several custom transformations which would be used solely by this
 * problem. Therefore we opt to implement it as a monolithic function instead.
 *
 * TODO: It would be nice to have a generic step ellipsoid function to
 * complement this one.
 */
#include <assert.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "coco_utilities.cpp"
#include "suite_bbob_legacy_code.cpp"

#include "transform_obj_norm_by_dim.cpp"
#include "transform_obj_penalize.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_blockrotation.cpp"
#include "transform_vars_permutation.cpp"
#include "transform_functions.cpp"

/**
 * @brief The BBOB step ellipsoid function.
 */
class f_step_ellipsoid_t final : public coco_problem_s {
public:
  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double raw(const double *x);

  /** @brief Creates the step ellipsoid function of the given instance (random seed). */
  f_step_ellipsoid_t(const size_t dimension, const long rseed, const coco_bbob_variant_t variant, const double fopt,
                     const double penalty_scale)
      : coco_problem_s(dimension, 1, 0), fopt(fopt), penalty_scale(penalty_scale) {
    static const double condition = 100;
    init_from_scalars("step ellipsoid function", -5.0, 5.0, 0);
    /* Allocate temporary storage */
    z_hat.resize(dimension);
    z.resize(dimension);
    xopt.resize(dimension);
    c1.resize(dimension);
    scales.resize(dimension);
    for (size_t i = 0; i < dimension; ++i) {
      c1[i] = sqrt(pow(condition / 10., static_cast<double>(i) / static_cast<double>(dimension - 1)));
      scales[i] = pow(condition, static_cast<double>(i) / (static_cast<double>(dimension) - 1.0));
    }
    bbob_compute_xopt(xopt.data(), rseed, dimension, variant);
    rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
    rot2 = bbob2009_compute_rotation(rseed, dimension);

    /* Compute best solution
     *
     * OME: Dirty hack for now because I did not want to invert the
     * transformations to find the best_parameter :/
     */
    for (size_t i = 0; i < number_of_variables; i++) {
      best_parameter[i] = xopt[i];
    }
    best_value[0] = fopt;
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = raw(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

private:
  std::vector<double> z_hat, z;
  std::vector<double> xopt;
  std::vector<double> c1;     /**< @brief sqrt((condition / 10)^(i / (n - 1))), scales the rows of rot2 */
  std::vector<double> scales; /**< @brief condition^(i / (n - 1)) */
  double fopt{}, penalty_scale{};
  coco_matrix_t rot1, rot2;
};

/**
 * @brief Implements the step ellipsoid function without connections to any COCO
 * structures.
 */
double f_step_ellipsoid_t::raw(const double *x) {

  static const double alpha = 10.0;
  double penalty = 0.0, x1;

  assert(number_of_variables > 1);

  for (size_t i = 0; i < number_of_variables; ++i) {
    double tmp = fabs(x[i]) - 5.0;
    if (tmp > 0.0)
      penalty += tmp * tmp;
  }
  penalty = penalty * (penalty_scale);
  for (size_t i = 0; i < number_of_variables; ++i) {
    z_hat[i] = 0.0;
    for (size_t j = 0; j < number_of_variables; ++j) {
      z_hat[i] += c1[i] * rot2(i, j) * (x[j] - xopt[j]);
    }
  }
  x1 = z_hat[0];

  for (size_t i = 0; i < number_of_variables; ++i) {
    if (fabs(z_hat[i]) > 0.5) /* TODO: Documentation: no fabs() in documentation */
      z_hat[i] = coco_double_round(z_hat[i]);
    else
      z_hat[i] = coco_double_round(alpha * z_hat[i]) / alpha;
  }

  rot1.multiply(z_hat.data(), z.data());

  /* Computation core */
  double result = 0.0;
  for (size_t i = 0; i < number_of_variables; ++i)
    result += scales[i] * z[i] * z[i];
  result = 0.1 * coco_double_max(fabs(x1) * 1.0e-4, result) + penalty + fopt;

  return result;
}

/**
 * @brief Creates the BBOB step ellipsoid problem.
 *
 * @note There is no separate basic allocate function.
 */
static coco_problem_ptr_t
f_step_ellipsoid_bbob_problem_allocate(const size_t function, const size_t dimension, const size_t instance,
                                       const long rseed, const coco_bbob_variant_t variant, const double penalty_scale,
                                       const char *problem_id_template, const char *problem_name_template) {

  coco_problem_ptr_t problem = std::make_unique<f_step_ellipsoid_t>(
      dimension, rseed, variant, bbob2009_compute_fopt(function, instance), penalty_scale);
  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("2-moderate");

  return problem;
}

/**
 * @brief The step ellipsoid function f7 of the bbob-largescale suite, with permuted block rotations.
 *
 * f(x) = 0.1 max(|ẑ_1| n 1e-4, f_ellipsoid^100(Q round(ẑ))) / n + f_pen(x) + fopt with ẑ = Λ^10 R (x - xopt), the
 * rounding to steps of 1 or 0.1 (for |ẑ_i| <= 0.5), and the permuted block rotations R and Q.
 */
class bbob_step_ellipsoid_permblockdiag_t final : public coco_problem_s {
public:
  bbob_step_ellipsoid_permblockdiag_t(const size_t function, const size_t dimension, const size_t instance,
                                      const long rseed, const coco_bbob_variant_t variant,
                                      const char *problem_id_template, const char *problem_name_template)
      : coco_problem_s(dimension, 1, 0), shift(bbob_xopt(rseed, dimension, variant).data(), dimension),
        rotation1(dimension, rseed + 2000000, rseed, rseed + 3000000), conditioning(10.0, dimension),
        round_step(10.0, dimension), rotation2(dimension, rseed + 4000000, rseed + 1000000, rseed + 5000000),
        scales(dimension), fopt(bbob2009_compute_fopt(function, instance)), z1(dimension), z2(dimension) {
    static const double condition = 100;
    for (size_t i = 0; i < dimension; ++i)
      scales[i] = pow(condition, static_cast<double>(i) / (static_cast<double>(dimension) - 1.0));
    init_from_scalars("step ellipsoid function", -5.0, 5.0, 0.0);

    /* The optimum of the ellipsoid is the origin, map it back through the transformations; the rounding keeps it */
    rotation2.inverse(best_parameter.data(), z1.data());
    round_step.apply(z1.data(), z1.data());
    conditioning.inverse(z1.data(), z1.data());
    rotation1.inverse(z1.data(), z2.data());
    shift.inverse(z2.data(), best_parameter.data());
    best_value[0] = (*this)(best_parameter.data());

    set_id(problem_id_template, function, instance, dimension);
    set_name(problem_name_template, function, instance, dimension);
    set_type("2-moderate");
  }

  double operator()(const double *x) {
    const double penalty = coco_boundary_penalty(x, smallest_values_of_interest.data(),
                                                 largest_values_of_interest.data(), number_of_variables);
    shift.apply(x, z1.data());
    rotation1.apply(z1.data(), z2.data());
    conditioning.apply(z2.data(), z2.data());
    /* multiplied by the dimension to counter-balance the normalization */
    const double zhat_1 = fabs(z2[0]) * static_cast<double>(number_of_variables);
    round_step.apply(z2.data(), z2.data());
    rotation2.apply(z2.data(), z1.data());

    double y = 0.0;
    for (size_t i = 0; i < number_of_variables; ++i)
      y += scales[i] * z1[i] * z1[i];
    y = 0.1 * coco_double_max(zhat_1 * 1.0e-4, y);
    return y * coco_norm_by_dim_factor(number_of_variables) + penalty_factor * penalty + fopt;
  }

  void evaluate_function(const double *x, double *y) override {
    y[0] = (*this)(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

private:
  static constexpr double penalty_factor = 1.0;

  const vars_shift_t shift;
  const vars_permblockdiag_t rotation1;
  const vars_conditioning_t conditioning;
  const vars_round_step_t round_step;
  const vars_permblockdiag_t rotation2;
  std::vector<double> scales; /**< @brief condition^(i / (n - 1)) */
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};
