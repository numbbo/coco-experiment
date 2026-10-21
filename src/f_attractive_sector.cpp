/**
 * @file f_attractive_sector.cpp
 * @brief Implementation of the attractive sector function and problem.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "bbob_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

/**
 * @brief The attractive sector function.
 */
class f_attractive_sector_t final : public coco_problem_s {
public:
  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double raw(std::span<const double> x) const;

  f_attractive_sector_t(const size_t number_of_variables, std::span<const double> xopt)
      : coco_problem_s(number_of_variables, 1, 0) {
    assert(xopt.size() == number_of_variables);
    init_from_scalars("attractive sector function", -5.0, 5.0, 0.0);
    set_id("%s_d%02zu", "attractive_sector", number_of_variables);

    this->xopt.assign(xopt.begin(), xopt.end());

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
  std::vector<double> xopt;
};

/**
 * @brief Implements the attractive sector function without connections to any
 * COCO structures.
 */
double f_attractive_sector_t::raw(std::span<const double> x) const {

  double result = 0.0;
  for (size_t i = 0; i < number_of_variables; ++i) {
    if (xopt[i] * x[i] > 0.0) {
      result += 100.0 * 100.0 * x[i] * x[i];
    } else {
      result += x[i] * x[i];
    }
  }
  return result;
}

/**
 * @brief The BBOB attractive sector function f6, also in the bbob-boxed suite.
 *
 * f(x) = T_osz(f_attractive_sector(Q Λ^10 R (x - xopt)))^0.9 + fopt.
 */
class bbob_attractive_sector_t final : public bbob_problem_t<bbob_attractive_sector_t> {
public:
  bbob_attractive_sector_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                           const coco_bbob_variant_t variant, const char *problem_id_template,
                           const char *problem_name_template)
      : bbob_attractive_sector_t(function, dimension, instance, rseed, bbob_xopt(rseed, dimension, variant),
                                 problem_id_template, problem_name_template) {}

  double operator()(std::span<const double> x) {
    shift.apply(x, z1);
    conditioned_rotation.apply(z1, z2);
    const double y = pow(coco_oscillate_objective(sector.raw(z2)), 0.9);
    return y + fopt;
  }

private:
  bbob_attractive_sector_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                           const std::vector<double> &xopt, const char *problem_id_template,
                           const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "2-moderate", problem_id_template, problem_name_template),
        shift(xopt),
        conditioned_rotation(conditioned_rotation_matrix(rseed, dimension), std::vector<double>(dimension)),
        sector(dimension, xopt), fopt(bbob2009_compute_fopt(function, instance)), z1(dimension), z2(dimension) {
    /* The optimum of the attractive sector function is the origin, map it back through the transformations */
    conditioned_rotation.inverse(best_parameter, z1);
    shift.inverse(z1, best_parameter);
    set_best_value_at_best_parameter();
  }

  /** @brief Returns Q Λ^10 R with the rotations Q and R. */
  static coco_matrix_t conditioned_rotation_matrix(const long rseed, const size_t dimension) {
    const coco_matrix_t rot1 = bbob2009_compute_rotation(rseed + 1000000, dimension);
    const coco_matrix_t rot2 = bbob2009_compute_rotation(rseed, dimension);
    std::vector<double> scales(dimension);
    for (size_t k = 0; k < dimension; ++k) {
      const double exponent = 1.0 * static_cast<int>(k) / (static_cast<double>(dimension) - 1.0);
      scales[k] = pow(sqrt(10.0), exponent);
    }
    return coco_matrix_product_with_diagonal(rot1, scales, rot2);
  }

  const vars_shift_t shift;
  const vars_affine_t conditioned_rotation;
  const f_attractive_sector_t sector;
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};

/**
 * @brief The attractive sector function f6 of the bbob-largescale suite, with permuted block rotations.
 *
 * f(x) = T_osz(f_attractive_sector(Q Λ^10 R (x - xopt)) / n)^0.9 + fopt with the permuted block rotations R and Q.
 */
class bbob_attractive_sector_permblockdiag_t final : public bbob_problem_t<bbob_attractive_sector_permblockdiag_t> {
public:
  bbob_attractive_sector_permblockdiag_t(const size_t function, const size_t dimension, const size_t instance,
                                         const long rseed, const coco_bbob_variant_t variant,
                                         const char *problem_id_template, const char *problem_name_template)
      : bbob_attractive_sector_permblockdiag_t(function, dimension, instance, rseed,
                                               bbob_xopt(rseed, dimension, variant), problem_id_template,
                                               problem_name_template) {}

  double operator()(std::span<const double> x) {
    shift.apply(x, z1);
    rotation1.apply(z1, z2);
    conditioning.apply(z2, z2);
    rotation2.apply(z2, z1);
    double y = sector.raw(z1) * coco_norm_by_dim_factor(number_of_variables);
    y = pow(coco_oscillate_objective(y), 0.9);
    return y + fopt;
  }

private:
  bbob_attractive_sector_permblockdiag_t(const size_t function, const size_t dimension, const size_t instance,
                                         const long rseed, const std::vector<double> &xopt,
                                         const char *problem_id_template, const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "2-moderate", problem_id_template, problem_name_template),
        shift(xopt), rotation1(dimension, rseed + 2000000, rseed, rseed + 3000000), conditioning(10.0, dimension),
        rotation2(dimension, rseed + 4000000, rseed + 1000000, rseed + 5000000), sector(dimension, xopt),
        fopt(bbob2009_compute_fopt(function, instance)), z1(dimension), z2(dimension) {
    /* The optimum of the attractive sector function is the origin, map it back through the transformations */
    rotation2.inverse(best_parameter, z1);
    conditioning.inverse(z1, z1);
    rotation1.inverse(z1, z2);
    shift.inverse(z2, best_parameter);
    set_best_value_at_best_parameter();
  }

  const vars_shift_t shift;
  const vars_permblockdiag_t rotation1;
  const vars_conditioning_t conditioning;
  const vars_permblockdiag_t rotation2;
  const f_attractive_sector_t sector;
  const double fopt;
  std::vector<double> z1, z2; /**< @brief The transformed decision values */
};
