/**
 * @file f_gallagher.cpp
 * @brief Implementation of the Gallagher function and problem.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "coco_utilities.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_obj_norm_by_dim.cpp"
#include "transform_obj_shift.cpp"
#include "transform_vars_blockrotation.cpp"
#include "transform_vars_permutation.cpp"
#include "transform_vars_permutation_helpers.cpp"
#include "transform_vars_scale.cpp"

#include "transform_vars_conditioning.cpp"
#include "transform_obj_oscillate.cpp"
#include "transform_obj_power.cpp"
#include "transform_obj_penalize.cpp"

#include "f_sphere.cpp"

/**
 * @brief Data of the large-scale Gallagher function: its sub-problems, one per peak.
 */
struct f_gallagher_shared_data_t {
  size_t number_of_peaks = 0;
  std::vector<coco_problem_ptr_t> sub_problems; /**< @brief Owned by the large-scale Gallagher problem */
};
/**
 * @brief A random permutation type for the Gallagher problem.
 *
 * Needed to create a random permutation that is compatible with the old code.
 */
struct f_gallagher_permutation_t {
  double value;
  size_t index;
};

/**
 * Comparison function used for sorting.
 */
static int f_gallagher_compare_doubles(const void *a, const void *b) {
  double temp = static_cast<const f_gallagher_permutation_t *>(a)->value -
                static_cast<const f_gallagher_permutation_t *>(b)->value;
  if (temp > 0)
    return 1;
  else if (temp < 0)
    return -1;
  else
    return 0;
}

/**
 * @brief The BBOB Gallagher function.
 */
class f_gallagher_t final : public coco_problem_s {
public:
  /** @brief Evaluates the function at x without any COCO bookkeeping. */
  double raw(const double *x);

  /**
   * @brief Creates the Gallagher function of the given instance (random seed) with its peaks, which lie in
   * [-5, 5] for the boxed variant and in [-4, 4] otherwise.
   */
  f_gallagher_t(const size_t dimension, const size_t number_of_peaks, const double penalty_scale, const long rseed,
                const coco_bbob_variant_t variant)
      : coco_problem_s(dimension, 1, 0) {
    /* problem_name will be overwritten below */
    init_from_scalars("Gallagher function", -5.0, 5.0, 0.0);

    /* Allocate temporary storage and space for the peaks */
    this->number_of_peaks = number_of_peaks;
    this->penalty_scale = penalty_scale;
    x_local = coco_matrix_t(number_of_peaks, dimension);
    arr_scales = coco_matrix_t(number_of_peaks, dimension);
    peak_values.resize(number_of_peaks);
    tmx.resize(dimension);

    const size_t peaks_21 = 21;
    const size_t peaks_101 = 101;

    double maxcondition = 1000.;
    /* maxcondition1 satisfies the old code and the doc but seems wrong in that it
     * is, with very high probability, not the largest condition level!!! */
    double maxcondition1 = 1000.;
    double fitvalues[2] = {1.1, 9.1};
    /* Parameters for generating local optima. In the old code, they are different in f21 and f22 */
    double b = 0, c = 0;
    /* Random permutation */
    std::vector<f_gallagher_permutation_t> rperm;

    if (number_of_peaks == peaks_101) {
      maxcondition1 = sqrt(maxcondition1);
      b = 10.;
      c = 5.;
    } else if (number_of_peaks == peaks_21) {
      b = 9.8;
      c = 4.9;
    } else {
      coco_error("f_gallagher_t(): '%zu' is a non-supported number of peaks", number_of_peaks);
    }
    rotation = bbob2009_compute_rotation(rseed, dimension);

    /* Initialize all the data of the inner problem */
    std::vector<double> random_numbers(number_of_peaks * dimension); /* This is large enough for all cases below */
    bbob2009_unif(random_numbers.data(), number_of_peaks - 1, rseed);
    rperm.resize(number_of_peaks - 1);
    for (size_t i = 0; i < number_of_peaks - 1; ++i) {
      rperm[i].value = random_numbers[i];
      rperm[i].index = i;
    }
    qsort(rperm.data(), number_of_peaks - 1, sizeof(rperm[0]), f_gallagher_compare_doubles);

    /* Random permutation */
    std::vector<double> arrCondition(number_of_peaks);
    arrCondition[0] = maxcondition1;
    peak_values[0] = 10;
    for (size_t i = 1; i < number_of_peaks; ++i) {
      arrCondition[i] =
          pow(maxcondition, static_cast<double>(rperm[i - 1].index) / (static_cast<double>(number_of_peaks - 2)));
      peak_values[i] =
          static_cast<double>(i - 1) / static_cast<double>(number_of_peaks - 2) * (fitvalues[1] - fitvalues[0]) +
          fitvalues[0];
    }

    rperm.resize(dimension);
    for (size_t i = 0; i < number_of_peaks; ++i) {
      bbob2009_unif(random_numbers.data(), dimension, rseed + static_cast<long>(1000 * i));
      for (size_t j = 0; j < dimension; ++j) {
        rperm[j].value = random_numbers[j];
        rperm[j].index = j;
      }
      qsort(rperm.data(), dimension, sizeof(rperm[0]), f_gallagher_compare_doubles);
      for (size_t j = 0; j < dimension; ++j) {
        arr_scales(i, j) = pow(arrCondition[i], /* Lambda^alpha_i from the doc */
                               (static_cast<double>(rperm[j].index)) / (static_cast<double>(dimension - 1)) - 0.5);
      }
    }

    bbob2009_unif(random_numbers.data(), dimension * number_of_peaks, rseed);

    /* SBOX-COST wants optima in [-5, 5], the other suites in [-4, 4]. */
    double scale_factor = 0.8;
    if (variant == coco_bbob_variant_t::boxed) {
      scale_factor = 1.0;
    }

    for (size_t i = 0; i < dimension; ++i) {
      best_parameter[i] = scale_factor * (b * random_numbers[i] - c);
      for (size_t j = 0; j < number_of_peaks; ++j) {
        x_local(j, i) = 0.;
        for (size_t k = 0; k < dimension; ++k) {
          x_local(j, i) += rotation(i, k) * (b * random_numbers[j * dimension + k] - c);
        }
        if (j == 0) {
          x_local(j, i) *= scale_factor;
        }
      }
    }

    /* Compute best solution */
    evaluate_function(best_parameter.data(), best_value.data());
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
  coco_matrix_t rotation;
  coco_matrix_t x_local; /**< @brief The location of peak i in row i (number_of_peaks x dimension) */
  coco_matrix_t arr_scales;
  size_t number_of_peaks{};
  double penalty_scale{};
  std::vector<double> peak_values;
  std::vector<double> tmx; /**< @brief The rotated decision values */
};

/**
 * @brief Implements the Gallagher function without connections to any COCO
 * structures.
 */
double f_gallagher_t::raw(const double *x) {
  double a = 0.1;
  double tmp2, f = 0., f_add, tmp, f_pen = 0., f_true = 0.;

  double fac = -0.5 / static_cast<double>(number_of_variables);

  /* Boundary handling */
  for (size_t i = 0; i < number_of_variables; ++i) {
    tmp = fabs(x[i]) - 5.;
    if (tmp > 0.) {
      f_pen += tmp * tmp;
    }
  }
  f_pen = f_pen * penalty_scale;
  f_add = f_pen;
  /* Transformation in search space */
  /* TODO: this should rather be done in f_gallagher */
  rotation.multiply(x, tmx.data());
  /* Computation core*/
  for (size_t i = 0; i < number_of_peaks; ++i) {
    const double *peak = x_local.row(i);
    const double *scales = arr_scales.row(i);
    tmp2 = 0.;
    for (size_t j = 0; j < number_of_variables; ++j) {
      tmp = (tmx[j] - peak[j]);
      tmp2 += scales[j] * tmp * tmp;
    }
    tmp2 = peak_values[i] * exp(fac * tmp2);
    f = coco_double_max(f, tmp2);
  }

  f = 10. - f;
  if (f > 0) {
    f_true = log(f) / a;
    f_true = pow(exp(f_true + 0.49 * (sin(f_true) + sin(0.79 * f_true))), a);
  } else if (f < 0) {
    f_true = log(-f) / a;
    f_true =
        -pow(exp(f_true + 0.49 * (sin(0.55 * f_true) + sin(0.31 * f_true))), a);
  } else
    f_true = f;

  f_true *= f_true;
  f_true += f_add;
  double result = f_true;
  return result;
}

/**
 * @brief Creates the BBOB Gallagher problem.
 *
 * @note There is no separate basic allocate function.
 */
static coco_problem_ptr_t f_gallagher_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                            const size_t instance, const long rseed,
                                                            const coco_bbob_variant_t variant,
                                                            const size_t number_of_peaks, const double penalty_scale,
                                                            const char *problem_id_template,
                                                            const char *problem_name_template) {

  coco_problem_ptr_t problem =
      std::make_unique<f_gallagher_t>(dimension, number_of_peaks, penalty_scale, rseed, variant);

  double fopt = bbob2009_compute_fopt(function, instance);
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("5-weakly-structured");

  return problem;
}

/* TODO: Functions used in/related to the large scale suite are below. Eventually either merge with above after the standard version is updated with the new approach or put what's below in a separate file */



/**
 * @brief A sub problem (one peak) of the large-scale Gallagher function.
 */
class f_gallagher_sub_t final : public coco_problem_s {
public:
  explicit f_gallagher_sub_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    init_from_scalars("gallagher_sub function", -5.0, 5.0, 0.0);
    set_id("%s_d%04zu", "gallagher_sub", number_of_variables);
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
 * @brief Allocates the basic gallagher sub problem.
 */
static coco_problem_ptr_t f_gallagher_sub_problem_allocate(const size_t number_of_variables) {
  return std::make_unique<f_gallagher_sub_t>(number_of_variables);
}

/**
 * @brief Implements the gallagher function without connections to any COCO
 * structures. Wassim: core to not conflict with raw for now
 */
static double f_gallagher_core(const double *x, size_t number_of_variables,
                               f_gallagher_shared_data_t *shared) {

  coco_problem_t *problem_i;
  double result = 0;
  double y, w_i;
  double maxf = DBL_MAX;

  for (size_t i = 0; i < shared->number_of_peaks; i++) {
    problem_i = shared->sub_problems[i].get();
    problem_i->evaluate_function(x, &y);
    if (i == 0) {
      w_i = 10;
    } else {
      w_i = 1.1 + 8.0 * ((static_cast<double>(i) + 1) - 2.0) / ((static_cast<double>(shared->number_of_peaks)) - 2.0);
    }
    y = w_i * exp(-1.0 / (2.0 * (static_cast<double>(number_of_variables))) *
                  y); /* Wassim: problem_i->evaluate_function is the sphere on
            a transformed coordiante system with conditioning */
    if (i == 0 || maxf < y) {
      maxf = y;
    }
  }
  result = 10.0 - maxf;
  return result;
}

/**
 * @brief The Gallagher function of the large-scale suite.
 *
 * Owns its sub problems (one per peak), which are stored in its shared data.
 */
class f_gallagher_ls_t final : public coco_problem_s {
public:
  f_gallagher_ls_t(const size_t number_of_variables, size_t number_of_peaks)
      : coco_problem_s(number_of_variables, 1, 0) {

    init_from_scalars("gallagher function", -5.0, 5.0, 0.0);
    shared.number_of_peaks = number_of_peaks;
    shared.sub_problems.resize(number_of_peaks);
    for (size_t peak_index = 0; peak_index < number_of_peaks; peak_index++) {
      shared.sub_problems[peak_index] = f_gallagher_sub_problem_allocate(number_of_variables);
    }
    set_id("%s_d%04zu", "gallagher", number_of_variables);
    best_value[0] = 0;
  }

  /**
   * @brief Uses the core function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_gallagher_core(x, number_of_variables, &shared);
    if (!(y[0] + 1e-13 >= best_value[0])) {
      coco_warning("x[0]= %f: %f < %f", x[0], y[0] + 1e-13, best_value[0]);
    }
    assert(y[0] + 1e-13 >= best_value[0]);
  }

  /** @brief The data shared with the transformations and sub problems of this function. */
  f_gallagher_shared_data_t &shared_data() { return shared; }

private:
  f_gallagher_shared_data_t shared;
};

/**
 * @brief Allocates the basic gallagher problem.
 */
static std::unique_ptr<f_gallagher_ls_t> f_gallagher_problem_allocate(const size_t number_of_variables,
                                                                      size_t number_of_peaks) {
  return std::make_unique<f_gallagher_ls_t>(number_of_variables, number_of_peaks);
}

static coco_problem_ptr_t f_gallagher_permblockdiag_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                          const size_t instance, const long rseed,
                                                                          size_t number_of_peaks,
                                                                          const char *problem_id_template,
                                                                          const char *problem_name_template) {

  double penalty_factor = 1.0;
  coco_problem_ptr_t *problem_i;
  f_gallagher_shared_data_t *shared;
  size_t nb_blocks;
  const size_t peaks_21 = 21;
  const size_t peaks_101 = 101;
  double a = 0.8, b = 0, c = 0;
  double first_condition = 0;
  double alpha_i;
  /* first_condition satisfies the old code and the doc but seems wrong in that
   * it is, with very high probability, not the largest condition level!!! */

  double fopt = bbob2009_compute_fopt(function, instance);
  if (number_of_peaks == peaks_101) {
    first_condition = 1000;
    b = 10.0;
    c = 5.0;
  } else if (number_of_peaks == peaks_21) {
    first_condition = 1000 * 1000;
    b = 9.8;
    c = 4.9;
  } else {
    coco_error("f_gallagher_permblockdiag_bbob_problem_allocate(): '%zu' is a "
               "non-supported number of peaks",
               number_of_peaks);
  }

  std::vector<size_t> block_sizes = coco_get_block_sizes(&nb_blocks, dimension, "bbob-largescale");


  const coco_block_diagonal_matrix_t B = coco_compute_blockrotation(rseed, block_sizes.data(), nb_blocks);

  std::unique_ptr<f_gallagher_ls_t> gallagher = f_gallagher_problem_allocate(dimension, number_of_peaks);
  shared = &gallagher->shared_data(); /* shortcut */
  coco_problem_ptr_t problem = std::move(gallagher);

  std::vector<double> alpha_i_vals(number_of_peaks - 1);
  for (size_t i = 0; i < number_of_peaks - 1; i++) {
    alpha_i_vals[i] = pow(1000, 2 * static_cast<double>(i) / (static_cast<double>(number_of_peaks - 2)));
  }
  std::vector<double> tmp_uniform(dimension * number_of_peaks);
  bbob2009_unif(tmp_uniform.data(), number_of_peaks - 1, rseed);
  /* random permutation of the alpha_i's to allow sampling without replacement*/
  std::vector<size_t> P_alpha_i(number_of_peaks - 1);
  coco_compute_permutation_from_sequence(P_alpha_i.data(), tmp_uniform.data(), number_of_peaks - 1);
  bbob2009_unif(tmp_uniform.data(), dimension * number_of_peaks, rseed);
  std::vector<double> best_param_before_rotation(dimension);
  std::vector<double> best_param_after_rotation(dimension);
  std::vector<double> tmp_uniform2(dimension);
  for (size_t peak_index = 0; peak_index < number_of_peaks; peak_index++) {
    problem_i = &(shared->sub_problems[peak_index]);

    /* compute transformation parameters: */
    /* y_i and block-rotate it once and for all */
    std::vector<double> y_i(dimension);
    for (size_t i = 0; i < dimension; i++) {
      y_i[i] = b * tmp_uniform[i + peak_index * dimension] - c;
      if (peak_index == 0) {
        y_i[i] *= a;
        best_param_before_rotation[i] = y_i[i];
      }
    }
    const std::vector<double> y_i_unrotated = y_i;
    B.multiply(y_i_unrotated.data(), y_i.data());
    if (peak_index == 0) {
      for (size_t i = 0; i < dimension; i++) {
        best_param_after_rotation[i] = y_i[i];
      }
    }

    if (peak_index == 0) {
      alpha_i = first_condition;
    } else {
      alpha_i = alpha_i_vals[P_alpha_i[peak_index - 1]]; /*already square-rooted */
    }

    /* P_Lambda the permutation */
    std::vector<size_t> P_Lambda(dimension); /* random permutation of the values in C_i */
    bbob2009_unif(tmp_uniform2.data(), dimension, rseed + static_cast<long>(1000 * peak_index));
    coco_compute_permutation_from_sequence(P_Lambda.data(), tmp_uniform2.data(), dimension);

    /* apply var transformations to sub problem*/
    *problem_i = transform_vars_scale(std::move(*problem_i), 1. / sqrt(sqrt(sqrt(alpha_i)))); /* sqrt( alpha^1/4) */
    *problem_i = transform_vars_conditioning(std::move(*problem_i), sqrt(alpha_i));
    /**problem_i = transform_vars_blockrotation(std::move(*problem_i), B_const, dimension,
     * block_sizes, nb_blocks);*/
    *problem_i = transform_vars_inverse_permutation(std::move(*problem_i), P_Lambda.data(), dimension);
    *problem_i = transform_vars_shift(std::move(*problem_i), y_i.data(), 0);
  }

  problem->evaluate_best_value(best_param_after_rotation.data());
  problem = transform_vars_blockrotation(std::move(problem), B);
  problem->set_best_parameter(best_param_before_rotation.data());

  /*transform global objective function*/
  problem = transform_obj_oscillate(std::move(problem));
  problem = transform_obj_power(std::move(problem), 2.0);
  problem = transform_obj_penalize(std::move(problem), penalty_factor);
  problem = transform_obj_shift(std::move(problem), fopt);

  problem->set_id(problem_id_template, function, instance, dimension);
  problem->set_name(problem_name_template, function, instance, dimension);
  problem->set_type("5-weakly-structured");

  return problem;
}
