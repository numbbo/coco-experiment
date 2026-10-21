/**
 * @file f_rastrigin.cpp
 * @brief Implementation of the Rastrigin function and problem.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "coco_utilities.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_norm_by_dim.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_affine.cpp"
#include "transform_vars_asymmetric.cpp"
#include "transform_vars_blockrotation.cpp"
#include "transform_vars_conditioning.cpp"
#include "transform_vars_oscillate.cpp"
#include "transform_vars_permutation.cpp"
#include "transform_vars_shift.cpp"
#include "transform_functions.cpp"

/**
 * @brief Implements the Rastrigin function without connections to any COCO
 * structures.
 */
static double f_rastrigin_raw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  size_t i = 0;
  double sum1 = 0.0, sum2 = 0.0;

  for (i = 0; i < number_of_variables; ++i) {
    sum1 += cos(coco_two_pi * x[i]);
    sum2 += x[i] * x[i];
  }
  if (coco_is_inf(sum2)) /* cos(inf) -> nan */
    return sum2;
  double result = 10.0 * (static_cast<double>(number_of_variables) - sum1) + sum2;

  return result;
}

/**
 * @brief The Rastrigin function.
 */
class f_rastrigin_t final : public coco_problem_s {
public:
  explicit f_rastrigin_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("Rastrigin function", -5.0, 5.0, 0.0);
    /* TODO: make sure the gradient is computed correctly for the rotated Rastrigin */
    set_id("%s_d%02zu", "rastrigin", number_of_variables);

    /* Compute best solution */
    evaluate_function(best_parameter, best_value);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    assert(number_of_objectives == 1);
    y[0] = f_rastrigin_raw(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /**
   * @brief Evaluates the gradient of the raw Rastrigin function.
   */
  void evaluate_gradient(std::span<const double> x, std::span<double> y) override {

    for (size_t i = 0; i < number_of_variables; ++i) {
      y[i] = 2.0 * (10. * coco_pi * sin(coco_two_pi * x[i]) + x[i]);
    }
  }
};

/**
 * @brief Allocates the basic Rastrigin problem.
 */
static coco_problem_ptr_t f_rastrigin_allocate(const size_t number_of_variables) {
  return std::make_unique<f_rastrigin_t>(number_of_variables);
}

/**
 * @brief The BBOB Rastrigin function f3, also in the bbob-boxed and bbob-largescale suites.
 *
 * f(x) = f_rastrigin(Λ^10 T_asy^0.2(T_osz(x - xopt))) + fopt, normalized by the dimension before fopt is added in
 * the large-scale suite.
 */
class bbob_rastrigin_t final : public coco_problem_s {
public:
  bbob_rastrigin_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                   const coco_bbob_variant_t variant, const char *problem_id_template,
                   const char *problem_name_template)
      : coco_problem_s(dimension, 1, 0), shift(bbob_xopt(rseed, dimension, variant)), asymmetric(0.2, dimension),
        conditioning(10.0, dimension), fopt(bbob2009_compute_fopt(function, instance)),
        normalize(variant == coco_bbob_variant_t::largescale), z(dimension) {
    init_from_scalars("Rastrigin function", -5.0, 5.0, 0.0);

    /* The optimum of the Rastrigin function is the origin, map it back through the transformations */
    conditioning.inverse(best_parameter, z);
    asymmetric.inverse(z, z);
    oscillate.inverse(z, z);
    shift.inverse(z, best_parameter);
    best_value[0] = (*this)(best_parameter);

    set_id(problem_id_template, function, instance, dimension);
    set_name(problem_name_template, function, instance, dimension);
    set_type("1-separable");
  }

  double operator()(std::span<const double> x) {
    shift.apply(x, z);
    oscillate.apply(z, z);
    asymmetric.apply(z, z);
    conditioning.apply(z, z);
    double y = f_rastrigin_raw(z);
    if (normalize)
      y *= coco_norm_by_dim_factor(number_of_variables);
    return y + fopt;
  }

  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    y[0] = (*this)(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

private:
  const vars_shift_t shift;
  const vars_oscillate_t oscillate{};
  const vars_asymmetric_t asymmetric;
  const vars_conditioning_t conditioning;
  const double fopt;
  const bool normalize;
  std::vector<double> z; /**< @brief The transformed decision values */
};

/**
 * @brief The BBOB rotated Rastrigin function f15, also in the bbob-boxed suite.
 *
 * f(x) = f_rastrigin(R Λ^10 Q T_asy^0.2(T_osz(R (x - xopt)))) + fopt.
 */
class bbob_rastrigin_rotated_t final : public coco_problem_s {
public:
  bbob_rastrigin_rotated_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                           const coco_bbob_variant_t variant, const char *problem_id_template,
                           const char *problem_name_template)
      : coco_problem_s(dimension, 1, 0), shift(bbob_xopt(rseed, dimension, variant)),
        rotation(bbob2009_compute_rotation(rseed + 1000000, dimension), std::vector<double>(dimension)),
        asymmetric(0.2, dimension),
        conditioned_rotation(conditioned_rotation_matrix(rseed, dimension), std::vector<double>(dimension)),
        fopt(bbob2009_compute_fopt(function, instance)), z1(dimension), z2(dimension) {
    init_from_scalars("Rastrigin function", -5.0, 5.0, 0.0);

    /* The optimum of the Rastrigin function is the origin, map it back through the transformations */
    conditioned_rotation.inverse(best_parameter, z1);
    asymmetric.inverse(z1, z1);
    oscillate.inverse(z1, z1);
    rotation.inverse(z1, z2);
    shift.inverse(z2, best_parameter);
    best_value[0] = (*this)(best_parameter);

    set_id(problem_id_template, function, instance, dimension);
    set_name(problem_name_template, function, instance, dimension);
    set_type("4-multi-modal");
  }

  double operator()(std::span<const double> x) {
    shift.apply(x, z1);
    rotation.apply(z1, z2);
    oscillate.apply(z2, z2);
    asymmetric.apply(z2, z2);
    conditioned_rotation.apply(z2, z1);
    return f_rastrigin_raw(z1) + fopt;
  }

  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    y[0] = (*this)(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

private:
  /** @brief Returns R Λ^10 Q with the rotations R and Q. */
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
  const vars_affine_t rotation;
  const vars_oscillate_t oscillate{};
  const vars_asymmetric_t asymmetric;
  const vars_affine_t conditioned_rotation;
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};

/**
 * @brief The rotated Rastrigin function f15 of the bbob-largescale suite, with permuted block rotations.
 *
 * f(x) = f_rastrigin(R Λ^10 Q T_asy^0.2(T_osz(R (x - xopt)))) / n + fopt with the permuted block rotations R and Q.
 */
class bbob_rastrigin_rotated_permblockdiag_t final : public coco_problem_s {
public:
  bbob_rastrigin_rotated_permblockdiag_t(const size_t function, const size_t dimension, const size_t instance,
                                         const long rseed, const coco_bbob_variant_t variant,
                                         const char *problem_id_template, const char *problem_name_template)
      : coco_problem_s(dimension, 1, 0), shift(bbob_xopt(rseed, dimension, variant)),
        rotation1(dimension, rseed + 3000000, rseed + 1000000, rseed + 4000000), asymmetric(0.2, dimension),
        rotation2(dimension, rseed + 5000000, rseed, rseed + 6000000), conditioning(10.0, dimension),
        fopt(bbob2009_compute_fopt(function, instance)), z1(dimension), z2(dimension) {
    init_from_scalars("Rastrigin function", -5.0, 5.0, 0.0);

    /* The optimum of the Rastrigin function is the origin, map it back through the transformations */
    rotation1.inverse(best_parameter, z1);
    conditioning.inverse(z1, z1);
    rotation2.inverse(z1, z2);
    asymmetric.inverse(z2, z2);
    oscillate.inverse(z2, z2);
    rotation1.inverse(z2, z1);
    shift.inverse(z1, best_parameter);
    best_value[0] = (*this)(best_parameter);

    set_id(problem_id_template, function, instance, dimension);
    set_name(problem_name_template, function, instance, dimension);
    set_type("4-multi-modal");
  }

  double operator()(std::span<const double> x) {
    shift.apply(x, z1);
    rotation1.apply(z1, z2);
    oscillate.apply(z2, z2);
    asymmetric.apply(z2, z2);
    rotation2.apply(z2, z1);
    conditioning.apply(z1, z1);
    rotation1.apply(z1, z2);
    const double y = f_rastrigin_raw(z2);
    return y * coco_norm_by_dim_factor(number_of_variables) + fopt;
  }

  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    y[0] = (*this)(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

private:
  const vars_shift_t shift;
  const vars_permblockdiag_t rotation1;
  const vars_oscillate_t oscillate{};
  const vars_asymmetric_t asymmetric;
  const vars_permblockdiag_t rotation2;
  const vars_conditioning_t conditioning;
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};

/**
 * @brief Computes xopt for constrained Rastrigin (alternative to
 * bbob2009_compute_xopt()) xopt is a vector of dim uniform random integers
 */
static void f_rastrigin_cons_compute_xopt(std::span<double> xopt, const long rseed) {
  bbob2009_unif(xopt, rseed);

  for (size_t i = 0; i < xopt.size(); ++i) {
    xopt[i] = 10 * xopt[i] - 5;
    xopt[i] = static_cast<int>(xopt[i]);
  }

  /* In case (0, ..., 0) is sampled, set xopt to a different value */
  if (coco_vector_is_zero(xopt))
    for (size_t i = 0; i < xopt.size(); ++i) {
      xopt[i] = static_cast<int>(i % 9) - 4;
    }
}

/**
 * @brief Creates the Rastrigin problem for the constrained BBOB suite.
 */
static coco_problem_ptr_t f_rastrigin_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                 const size_t instance, const long rseed,
                                                                 const char *problem_id_template,
                                                                 const char *problem_name_template) {

  coco_problem_ptr_t problem;

  std::vector<double> xopt(dimension);
  f_rastrigin_cons_compute_xopt(xopt, rseed);
  double fopt = bbob2009_compute_fopt(function, instance);

  problem = f_rastrigin_allocate(dimension);
  problem = transform_vars_shift(std::move(problem), xopt, 0);
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("1-separable");

  return problem;
}

/**
 * @brief Creates the rotated Rastrigin problem for the constrained BBOB suite.
 */
static coco_problem_ptr_t f_rastrigin_rotated_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                         const size_t instance, const long rseed,
                                                                         const char *problem_id_template,
                                                                         const char *problem_name_template) {

  coco_problem_ptr_t problem;

  coco_matrix_t M(dimension, dimension);
  std::vector<double> b(dimension);

  const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
  M = rot1;

  std::vector<double> xopt(dimension);
  f_rastrigin_cons_compute_xopt(xopt, rseed);
  double fopt = bbob2009_compute_fopt(function, instance);

  problem = f_rastrigin_allocate(dimension);
  problem = transform_vars_shift(std::move(problem), xopt, 0);
  problem = transform_vars_affine(std::move(problem), M, b); /* Rotate after shifting so R dot xopt is also rotated */
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("1-separable");

  return problem;
}
