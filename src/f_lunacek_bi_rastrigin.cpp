/**
 * @file f_lunacek_bi_rastrigin.cpp
 * @brief Implementation of the Lunacek bi-Rastrigin function and problem.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "f_sphere.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_norm_by_dim.cpp"
#include "transform_obj_scale.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_x_hat_generic.cpp"

/**
 * @brief The BBOB Lunacek bi-Rastrigin function.
 */
struct f_lunacek_bi_rastrigin_t final : public coco_problem_s {
  std::vector<double> x_hat, z;
  std::vector<double> xopt;
  double fopt{};
  coco_matrix_t rot1, rot2;
  long rseed{};

  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double raw(const double *x);

  explicit f_lunacek_bi_rastrigin_t(const size_t dimension) : coco_problem_s(dimension, 1, 0) {
    init_from_scalars("Lunacek's bi-Rastrigin function", -5.0, 5.0, 0.0);
    /* Allocate temporary storage and space for the rotation matrices */
    x_hat.resize(dimension);
    z.resize(dimension);
    xopt.resize(dimension);
    rot1 = coco_matrix_t(dimension, dimension);
    rot2 = coco_matrix_t(dimension, dimension);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = raw(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Implements the Lunacek bi-Rastrigin function without connections to
 * any COCO structures.
 */
double f_lunacek_bi_rastrigin_t::raw(const double *x) {
  double result;
  static const double condition = 100.;
  size_t i, j;
  double penalty = 0.0;
  static const double mu0 = 2.5;
  static const double d = 1.;
  const double s = 1. - 0.5 / (sqrt((double)(number_of_variables + 20)) - 4.1);
  const double mu1 = -sqrt((mu0 * mu0 - d) / s);
  double sum1 = 0., sum2 = 0., sum3 = 0.;

  assert(number_of_variables > 1);

  if (coco_vector_contains_nan(x, number_of_variables))
    return NAN;

  for (i = 0; i < number_of_variables; ++i) {
    double tmp;
    tmp = fabs(x[i]) - 5.0;
    if (tmp > 0.0)
      penalty += tmp * tmp;
  }

  /* x_hat */
  for (i = 0; i < number_of_variables; ++i) {
    x_hat[i] = 2. * x[i];
    if (xopt[i] < 0.) {
      x_hat[i] *= -1.;
    }
  }

  std::vector<double> tmpvect(number_of_variables);
  /* affine transformation */
  for (i = 0; i < number_of_variables; ++i) {
    double c1;
    tmpvect[i] = 0.0;
    c1 = pow(sqrt(condition), ((double)i) / (double)(number_of_variables - 1));
    for (j = 0; j < number_of_variables; ++j) {
      tmpvect[i] += c1 * rot2[i][j] * (x_hat[j] - mu0);
    }
  }
  for (i = 0; i < number_of_variables; ++i) {
    z[i] = 0;
    for (j = 0; j < number_of_variables; ++j) {
      z[i] += rot1[i][j] * tmpvect[j];
    }
  }
  /* Computation core */
  for (i = 0; i < number_of_variables; ++i) {
    sum1 += (x_hat[i] - mu0) * (x_hat[i] - mu0);
    sum2 += (x_hat[i] - mu1) * (x_hat[i] - mu1);
    sum3 += cos(2 * coco_pi * z[i]);
  }
  result = coco_double_min(sum1, d * (double)number_of_variables + s * sum2) +
           10. * ((double)number_of_variables - sum3) + 1e4 * penalty;

  return result;
}

/**
 * @brief Creates the BBOB Lunacek bi-Rastrigin problem.
 *
 * @note There is no separate basic allocate function.
 */
static coco_problem_ptr_t f_lunacek_bi_rastrigin_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                       const size_t instance, const long rseed,
                                                                       const char *problem_id_template,
                                                                       const char *problem_name_template) {

  std::unique_ptr<f_lunacek_bi_rastrigin_t> owned = std::make_unique<f_lunacek_bi_rastrigin_t>(dimension);
  f_lunacek_bi_rastrigin_t *lunacek = owned.get(); /* Stays valid, the transformations take ownership of it */
  coco_problem_ptr_t problem = std::move(owned);

  const double mu0 = 2.5;

  double fopt;
  size_t i;

  lunacek->rseed = rseed;

  lunacek->fopt = bbob2009_compute_fopt(24, instance);
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    sbox_cost_compute_xopt(lunacek->xopt.data(), rseed, dimension);
  } else {
    bbob2009_compute_xopt(lunacek->xopt.data(), rseed, dimension);
  }
  bbob2009_compute_rotation(lunacek->rot1.rows(), rseed + 1000000, dimension);
  bbob2009_compute_rotation(lunacek->rot2.rows(), rseed, dimension);

  /* Compute best solution */
  std::vector<double> tmpvect(dimension);
  bbob2009_gauss(tmpvect.data(), dimension, rseed);
  for (i = 0; i < dimension; ++i) {
    lunacek->xopt[i] = 0.5 * mu0;
    if (tmpvect[i] < 0.0) {
      lunacek->xopt[i] *= -1.0;
    }
    problem->best_parameter[i] = lunacek->xopt[i];
  }
  problem->evaluate_function(problem->best_parameter.data(), problem->best_value.data());

  fopt = bbob2009_compute_fopt(function, instance);
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("5-weakly-structured");

  return problem;
}

/* Functions used in/related to the large scale suite are below. Eventually either merge with above after the standard
 * version is updated with the new approach or put what's below in a separate file*/

/**
 * @brief A sub problem of the large-scale Lunacek bi-Rastrigin function (a sphere function).
 */
struct f_lunacek_bi_rastrigin_sub_t final : public coco_problem_s {
  explicit f_lunacek_bi_rastrigin_sub_t(const size_t number_of_variables)
      : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("lunacek_bi_rastrigin_sub function", -5.0, 5.0, 0.0);
    set_id("%s_d%04lu", "lunacek_bi_rastrigin_sub", number_of_variables);
    evaluate_function(best_parameter.data(), best_value.data());
  }

  /**
   * @brief Uses the core function to evaluate the sub problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_sphere_raw(x, number_of_variables);
  }
};

/**
 * @brief Allocates the basic lunacek_bi_rastrigin sub problem.
 */
static coco_problem_ptr_t f_lunacek_bi_rastrigin_sub_problem_allocate(const size_t number_of_variables) {
  return std::make_unique<f_lunacek_bi_rastrigin_sub_t>(number_of_variables);
}

/**
 * @brief Implements the lunacek_bi_rastrigin function without connections to
 * any COCO structures. Wassim: core to not conflict with raw for now
 */
static double
f_lunacek_bi_rastrigin_core(const double *x, const size_t number_of_variables,
                            f_lunacek_bi_rastrigin_shared_data_t *shared) {

  coco_problem_t *problem_sub_mu0, *problem_sub_mu1;
  size_t i;
  double result = 0.0;
  double y0, y1;
  double *x_hat;

  x_hat = shared->x_hat.data();

  problem_sub_mu0 = shared->sub_problem_mu0.get();
  problem_sub_mu0->evaluate_function(x_hat, &y0);

  problem_sub_mu1 = shared->sub_problem_mu1.get();
  problem_sub_mu1->evaluate_function(x_hat, &y1);

  result += (double)number_of_variables;

  for (i = 0; i < number_of_variables; i++) {
    result -= cos(2 * coco_pi * x[i]);
  }
  result *= 10.0;
  result += coco_double_min(y0, y1);

  return result;
}

/**
 * @brief The Lunacek bi-Rastrigin function of the large-scale suite.
 *
 * Shares x_hat with transform_vars_x_hat_generic, which sets it. Owns the two sub problems.
 */
struct f_lunacek_bi_rastrigin_ls_t final : public coco_problem_s {
  f_lunacek_bi_rastrigin_shared_data_t shared;

  explicit f_lunacek_bi_rastrigin_ls_t(const size_t number_of_variables)
      : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("lunacek_bi_rastrigin function", -5.0, 5.0, 0.0);
    shared.x_hat.resize(number_of_variables);
    set_id("%s_d%04lu", "lunacek_bi_rastrigin", number_of_variables);

    /* Compute the best solution later once the sub-problems are well defined */
    best_value[0] = 0; /* Manh: set default value to avoid assert() in transformation later*/
  }

  /**
   * @brief Uses the core function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_lunacek_bi_rastrigin_core(x, number_of_variables, &shared);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Allocates the basic lunacek_bi_rastrigin problem.
 */
static std::unique_ptr<f_lunacek_bi_rastrigin_ls_t>
f_lunacek_bi_rastrigin_problem_allocate(const size_t number_of_variables) {
  return std::make_unique<f_lunacek_bi_rastrigin_ls_t>(number_of_variables);
}

/**
 * @brief Creates the BBOB large scale suite Lunacek bi-Rastrigin problem.
 */
static coco_problem_ptr_t f_lunacek_bi_rastrigin_permblockdiag_bbob_problem_allocate(
    const size_t function, const size_t dimension, const size_t instance, const long rseed,
    const char *problem_id_template, const char *problem_name_template) {
  size_t i;
  double fopt;
  double penalty_factor = 1e4;
  coco_problem_ptr_t problem;
  coco_problem_ptr_t *sub_problem_tmp;

  size_t nb_blocks1, nb_blocks2;
  size_t swap_range1, swap_range2;
  size_t nb_swaps1, nb_swaps2;

  double condition = 100.0;
  double mu0, mu1, d = 1.0, s;

  fopt = bbob2009_compute_fopt(function, instance);
  s = 1. - 0.5 / (sqrt((double)(dimension + 20)) - 4.1);
  mu0 = 2.5;
  mu1 = -sqrt((mu0 * mu0 - d) / s);

  std::vector<size_t> block_sizes1 = coco_get_block_sizes(&nb_blocks1, dimension, "bbob-largescale");
  std::vector<size_t> block_sizes2 = coco_get_block_sizes(&nb_blocks2, dimension, "bbob-largescale");
  swap_range1 = coco_get_swap_range(dimension, "bbob-largescale");
  swap_range2 = coco_get_swap_range(dimension, "bbob-largescale");
  nb_swaps1 = coco_get_nb_swaps(dimension, "bbob-largescale");
  nb_swaps2 = coco_get_nb_swaps(dimension, "bbob-largescale");

  coco_matrix_t B1 = coco_block_matrix(dimension, block_sizes1.data(), nb_blocks1);
  coco_matrix_t B2 = coco_block_matrix(dimension, block_sizes2.data(), nb_blocks2);
  coco_compute_blockrotation(B1.rows(), rseed + 1000000, dimension, block_sizes1.data(),
                             nb_blocks1);
  coco_compute_blockrotation(B2.rows(), rseed, dimension, block_sizes2.data(), nb_blocks2);

  std::vector<size_t> P11(dimension);
  std::vector<size_t> P12(dimension);
  std::vector<size_t> P21(dimension);
  std::vector<size_t> P22(dimension);
  coco_compute_truncated_uniform_swap_permutation(
      P11.data(), rseed + 3000000, dimension, nb_swaps1, swap_range1);
  coco_compute_truncated_uniform_swap_permutation(
      P12.data(), rseed + 4000000, dimension, nb_swaps1, swap_range1);
  coco_compute_truncated_uniform_swap_permutation(
      P21.data(), rseed + 5000000, dimension, nb_swaps2, swap_range2);
  coco_compute_truncated_uniform_swap_permutation(
      P22.data(), rseed + 6000000, dimension, nb_swaps2, swap_range2);

  std::unique_ptr<f_lunacek_bi_rastrigin_ls_t> lunacek = f_lunacek_bi_rastrigin_problem_allocate(dimension);
  f_lunacek_bi_rastrigin_ls_t *raw_problem = lunacek.get();
  problem = std::move(lunacek);

  std::vector<double> mu0_vector(dimension);
  std::vector<double> mu1_vector(dimension);
  for (i = 0; i < dimension; i++) {
    mu0_vector[i] = mu0;
    mu1_vector[i] = mu1;
  }
  /* allocate sub-problems */
  raw_problem->shared.sub_problem_mu0 =
      f_lunacek_bi_rastrigin_sub_problem_allocate(dimension);
  raw_problem->shared.sub_problem_mu1 =
      f_lunacek_bi_rastrigin_sub_problem_allocate(dimension);

  /* set sign_vector */
  std::vector<double> tmp_normal(dimension);
  bbob2009_gauss(tmp_normal.data(), dimension, rseed);
  std::vector<double> sign_vector(dimension); /* Wassim sign vector designate the 1^+_- vector*/
  for (i = 0; i < dimension; i++) { /* set sign(x_opt)*/
    if (tmp_normal[i] < 0.0) {      /* Wassim: noraml is used here but unif for Schweffel!!! */
      sign_vector[i] = -1.0;
    } else {
      sign_vector[i] = 1.0;
    }
    coco_debug("An arbitrary message to suppress errors on Windows, see #1909. Please ignore accordingly.");
  }

  /* apply transformations to sub-problems */
  sub_problem_tmp = &raw_problem->shared.sub_problem_mu0;
  *sub_problem_tmp = transform_vars_shift(std::move(*sub_problem_tmp), mu0_vector.data(), 0);

  sub_problem_tmp = &raw_problem->shared.sub_problem_mu1;
  *sub_problem_tmp = transform_vars_shift(std::move(*sub_problem_tmp), mu1_vector.data(), 0);
  *sub_problem_tmp = transform_obj_scale(std::move(*sub_problem_tmp), s); /* Manh: need to use s */
  *sub_problem_tmp = transform_obj_shift(std::move(*sub_problem_tmp), d * (double)dimension);

  /* transformations on main problem */
  problem = transform_vars_permutation(std::move(problem), P22.data(), dimension);
  problem = transform_vars_blockrotation(std::move(problem), B1.rows(), dimension, block_sizes1.data(), nb_blocks1);
  problem = transform_vars_permutation(std::move(problem), P21.data(), dimension);
  problem = transform_vars_conditioning(std::move(problem), condition);

  problem = transform_vars_permutation(std::move(problem), P12.data(), dimension);
  problem = transform_vars_blockrotation(std::move(problem), B2.rows(), dimension, block_sizes2.data(), nb_blocks2);
  problem = transform_vars_permutation(std::move(problem), P11.data(), dimension);
  problem = transform_vars_shift(std::move(problem), mu0_vector.data(), 0);
  problem = transform_vars_x_hat_generic(std::move(problem), sign_vector.data(), &raw_problem->shared);

  problem = transform_obj_norm_by_dim(std::move(problem));
  problem = transform_obj_penalize(std::move(problem), penalty_factor);
  problem = transform_obj_shift(std::move(problem), fopt);

  /* set best_parameter and best value*/
  for (i = 0; i < dimension; i++) {                          /* Wassim: to silence warning about best_parameter*/
    problem->best_parameter[i] = 0.5 * mu0 * sign_vector[i]; /* TODO: Documentation no 0.5 in documentation! */
  }

  /* f_lunacek_bi_rastrigin_evaluate_core(problem.get(), problem->best_parameter, problem->best_value);
  printf("\n %f , x_opt[0]= %f\n", problem->best_value[0], problem->best_parameter[0]);*/ /* Wassim: for testing purposes, might end up being the one kept though*/

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("5-weakly-structured");

  return problem;
}
