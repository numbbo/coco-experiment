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
#include "transform_obj_norm_by_dim.cpp"

#include "transform_vars_conditioning.cpp"
#include "transform_obj_oscillate.cpp"
#include "transform_obj_power.cpp"
#include "transform_obj_penalize.cpp"

#include "f_sphere.cpp"
#include "transform_vars_gallagher_blockrotation.cpp"
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
 * @brief Data type for the Gallagher problem.
 */
struct f_gallagher_data_t {
  long rseed;
  std::vector<double> xopt;
  coco_matrix_t rotation, x_local, arr_scales;
  size_t number_of_peaks;
  double penalty_scale;
  std::vector<double> peak_values;
};

/**
 * Comparison function used for sorting.
 */
static int f_gallagher_compare_doubles(const void *a, const void *b) {
  double temp = (*(f_gallagher_permutation_t *)a).value - (*(f_gallagher_permutation_t *)b).value;
  if (temp > 0)
    return 1;
  else if (temp < 0)
    return -1;
  else
    return 0;
}

/**
 * @brief Implements the Gallagher function without connections to any COCO
 * structures.
 */
static double f_gallagher_raw(const double *x, const size_t number_of_variables,
                              f_gallagher_data_t *data) {
  size_t i, j; /* Loop over dim */
  double a = 0.1;
  double tmp2, f = 0., f_add, tmp, f_pen = 0., f_true = 0.;
  double fac;
  double result;

  if (coco_vector_contains_nan(x, number_of_variables))
    return NAN;

  fac = -0.5 / (double)number_of_variables;

  /* Boundary handling */
  for (i = 0; i < number_of_variables; ++i) {
    tmp = fabs(x[i]) - 5.;
    if (tmp > 0.) {
      f_pen += tmp * tmp;
    }
  }
  f_pen = f_pen * data->penalty_scale;
  f_add = f_pen;
  /* Transformation in search space */
  /* TODO: this should rather be done in f_gallagher */
  std::vector<double> tmx(number_of_variables);
  for (i = 0; i < number_of_variables; i++) {
    tmx[i] = 0;
    for (j = 0; j < number_of_variables; ++j) {
      tmx[i] += data->rotation[i][j] * x[j];
    }
  }
  /* Computation core*/
  for (i = 0; i < data->number_of_peaks; ++i) {
    tmp2 = 0.;
    for (j = 0; j < number_of_variables; ++j) {
      tmp = (tmx[j] - data->x_local[j][i]);
      tmp2 += data->arr_scales[i][j] * tmp * tmp;
    }
    tmp2 = data->peak_values[i] * exp(fac * tmp2);
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
  result = f_true;
  return result;
}

/**
 * @brief The BBOB Gallagher function.
 */
struct f_gallagher_t final : public coco_problem_s {
  f_gallagher_data_t params;

  f_gallagher_t(const size_t dimension, const size_t number_of_peaks, const double penalty_scale)
      : coco_problem_s(dimension, 1, 0) {
    /* problem_name and best_parameter will be overwritten below */
    coco_problem_init_from_scalars(this, "Gallagher function", -5.0, 5.0, 0.0);

    /* Allocate temporary storage and space for the rotation matrices */
    params.number_of_peaks = number_of_peaks;
    params.penalty_scale = penalty_scale;
    params.xopt.resize(dimension);
    params.rotation = coco_matrix_t(dimension, dimension);
    params.x_local = coco_matrix_t(dimension, number_of_peaks);
    params.arr_scales = coco_matrix_t(number_of_peaks, dimension);
    params.peak_values.resize(number_of_peaks);
  }

  /**
   * @brief Uses the raw function to evaluate the COCO problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_gallagher_raw(x, number_of_variables, &params);
    assert(y[0] + 1e-13 >= best_value[0]);
  }
};

/**
 * @brief Creates the BBOB Gallagher problem.
 *
 * @note There is no separate basic allocate function.
 */
static coco_problem_ptr_t f_gallagher_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                            const size_t instance, const long rseed, const void *args,
                                                            const char *problem_id_template,
                                                            const char *problem_name_template) {

  f_gallagher_data_t *data;
  coco_problem_ptr_t problem;

  const size_t peaks_21 = 21;
  const size_t peaks_101 = 101;

  double fopt;
  size_t i, j, k;
  double maxcondition = 1000.;
  /* maxcondition1 satisfies the old code and the doc but seems wrong in that it
   * is, with very high probability, not the largest condition level!!! */
  double maxcondition1 = 1000.;
  double fitvalues[2] = {1.1, 9.1};
  /* Parameters for generating local optima. In the old code, they are different in f21 and f22 */
  double b = 0, c = 0;
  /* Random permutation */
  std::vector<f_gallagher_permutation_t> rperm;

  f_gallagher_args_t *f_gallagher_args;
  f_gallagher_args = ((f_gallagher_args_t *)args);

  size_t number_of_peaks = f_gallagher_args->number_of_peaks;
  double penalty_scale = f_gallagher_args->penalty_scale;
  problem = std::make_unique<f_gallagher_t>(dimension, number_of_peaks, penalty_scale);
  data = &static_cast<f_gallagher_t *>(problem.get())->params;

  if (data->number_of_peaks == peaks_101) {
    maxcondition1 = sqrt(maxcondition1);
    b = 10.;
    c = 5.;
  } else if (data->number_of_peaks == peaks_21) {
    b = 9.8;
    c = 4.9;
  } else {
    coco_error("f_gallagher_bbob_problem_allocate(): '%lu' is a non-supported number of peaks",
               (unsigned long)data->number_of_peaks);
  }
  data->rseed = rseed;
  bbob2009_compute_rotation(data->rotation.rows(), rseed, dimension);

  /* Initialize all the data of the inner problem */
  std::vector<double> random_numbers(data->number_of_peaks * dimension); /* This is large enough for all cases below */
  bbob2009_unif(random_numbers.data(), data->number_of_peaks - 1, data->rseed);
  rperm.resize(data->number_of_peaks - 1);
  for (i = 0; i < data->number_of_peaks - 1; ++i) {
    rperm[i].value = random_numbers[i];
    rperm[i].index = i;
  }
  qsort(rperm.data(), data->number_of_peaks - 1, sizeof(rperm[0]), f_gallagher_compare_doubles);

  /* Random permutation */
  std::vector<double> arrCondition(data->number_of_peaks);
  arrCondition[0] = maxcondition1;
  data->peak_values[0] = 10;
  for (i = 1; i < data->number_of_peaks; ++i) {
    arrCondition[i] = pow(maxcondition, (double)(rperm[i - 1].index) / ((double)(data->number_of_peaks - 2)));
    data->peak_values[i] =
        (double)(i - 1) / (double)(data->number_of_peaks - 2) * (fitvalues[1] - fitvalues[0]) + fitvalues[0];
  }

  rperm.resize(dimension);
  for (i = 0; i < data->number_of_peaks; ++i) {
    bbob2009_unif(random_numbers.data(), dimension, data->rseed + (long)(1000 * i));
    for (j = 0; j < dimension; ++j) {
      rperm[j].value = random_numbers[j];
      rperm[j].index = j;
    }
    qsort(rperm.data(), dimension, sizeof(rperm[0]), f_gallagher_compare_doubles);
    for (j = 0; j < dimension; ++j) {
      data->arr_scales[i][j] = pow(arrCondition[i], /* Lambda^alpha_i from the doc */
                                   ((double)rperm[j].index) / ((double)(dimension - 1)) - 0.5);
    }
  }

  bbob2009_unif(random_numbers.data(), dimension * data->number_of_peaks, data->rseed);

  /* SBOX-COST wants optima in [-5, 5], the other suites in [-4, 4]. */
  double scale_factor = 0.8;
  if (coco_strfind(problem_name_template, "bbob-boxed suite problem") >= 0) {
    scale_factor = 1.0;
  }

  for (i = 0; i < dimension; ++i) {
    data->xopt[i] = scale_factor * (b * random_numbers[i] - c);
    problem->best_parameter[i] = scale_factor * (b * random_numbers[i] - c);
    for (j = 0; j < data->number_of_peaks; ++j) {
      data->x_local[i][j] = 0.;
      for (k = 0; k < dimension; ++k) {
        data->x_local[i][j] +=
            data->rotation[i][k] * (b * random_numbers[j * dimension + k] - c);
      }
      if (j == 0) {
        data->x_local[i][j] *= scale_factor;
      }
    }
  }

  /* Compute best solution */
  problem->evaluate_function(problem->best_parameter.data(), problem->best_value.data());

  fopt = bbob2009_compute_fopt(function, instance);
  problem = transform_obj_shift(std::move(problem), fopt);

  coco_problem_set_id(problem.get(), problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem.get(), problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem.get(), "5-weakly-structured");

  return problem;
}

/* TODO: Functions used in/related to the large scale suite are below. Eventually either merge with above after the standard version is updated with the new approach or put what's below in a separate file */



/**
 * @brief A sub problem (one peak) of the large-scale Gallagher function.
 */
struct f_gallagher_sub_t final : public coco_problem_s {
  explicit f_gallagher_sub_t(const size_t number_of_variables) : coco_problem_s(number_of_variables, 1, 0) {
    coco_problem_init_from_scalars(this, "gallagher_sub function", -5.0, 5.0, 0.0);
    coco_problem_set_id(this, "%s_d%04lu", "gallagher_sub", number_of_variables);

    /* Compute best solution */
    /*f_gallagher_sub_evaluate_core(problem, problem->best_parameter,
     * problem->best_value);*/
  }

  /**
   * @brief Uses the core function to evaluate the sub problem.
   */
  void evaluate_function(const double *x, double *y) override {
    assert(number_of_objectives == 1);
    y[0] = f_sphere_raw(x, number_of_variables);
    /*assert(y[0] + 1e-13 >= problem->best_value[0]);*/
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
  size_t i, j;
  double maxf = DBL_MAX;
  std::vector<double> x_local(number_of_variables);

  for (i = 0; i < shared->number_of_peaks; i++) {
    for (j = 0; j < number_of_variables; j++) {
      x_local[j] = x[j]; /*shared->rotated_x[j];*/
    }
    problem_i = shared->sub_problems[i].get();
    problem_i->evaluate_function(x_local.data(), &y);
    if (i == 0) {
      w_i = 10;
    } else {
      w_i = 1.1 + 8.0 * (((double)i + 1) - 2.0) / (((double)shared->number_of_peaks) - 2.0);
    }
    y = w_i *
        exp(-1.0 / (2.0 * ((double)number_of_variables)) * y); /* Wassim: problem_i->evaluate_function is the sphere on
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
 * Owns its sub problems (one per peak), which are stored in its shared data together with the block rotation.
 */
struct f_gallagher_ls_t final : public coco_problem_s {
  f_gallagher_shared_data_t shared;

  f_gallagher_ls_t(const size_t number_of_variables, size_t number_of_peaks)
      : coco_problem_s(number_of_variables, 1, 0) {
    size_t peak_index;

    coco_problem_init_from_scalars(this, "gallagher function", -5.0, 5.0, 0.0);
    shared.number_of_peaks = number_of_peaks;
    shared.sub_problems.resize(number_of_peaks);
    for (peak_index = 0; peak_index < number_of_peaks; peak_index++) {
      shared.sub_problems[peak_index] = f_gallagher_sub_problem_allocate(number_of_variables);
    }
    shared.rotated_x.resize(number_of_variables);
    shared.block_size_map.resize(number_of_variables);
    shared.first_non_zero_map.resize(number_of_variables);
    shared.block_sizes.resize(number_of_variables);
    coco_problem_set_id(this, "%s_d%04lu", "gallagher", number_of_variables);
    best_value[0] = 0;
    /* Compute best solution */
    /*f_gallagher_evaluate_core(problem, problem->best_parameter,
     * problem->best_value);*/
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

  size_t i;
  size_t peak_index;
  double fopt;
  double penalty_factor = 1.0;
  coco_problem_ptr_t *problem_i;
  coco_problem_ptr_t rotation_problem;
  f_gallagher_shared_data_t *shared;
  /*size_t *P1, *P2;*/
  size_t nb_blocks;
  size_t idx_blocksize, current_blocksize,
      next_bs_change; /* needed for the rotated y_i*/
  /*size_t swap_range;
  size_t nb_swaps;*/
  const size_t peaks_21 = 21;
  const size_t peaks_101 = 101;
  double a = 0.8, b = 0, c = 0;
  double first_condition = 0;
  double alpha_i;
  /* first_condition satisfies the old code and the doc but seems wrong in that
   * it is, with very high probability, not the largest condition level!!! */

  fopt = bbob2009_compute_fopt(function, instance);
  if (number_of_peaks == peaks_101) {
    first_condition = 1000;
    b = 10.0;
    c = 5.0;
  } else if (number_of_peaks == peaks_21) {
    first_condition = 1000 * 1000;
    b = 9.8;
    c = 4.9;
  } else {
    coco_error("f_gallagher_permblockdiag_bbob_problem_allocate(): '%lu' is a "
               "non-supported number of peaks",
               number_of_peaks);
  }

  std::vector<size_t> block_sizes = coco_get_block_sizes(&nb_blocks, dimension, "bbob-largescale");

  /*swap_range = coco_get_swap_range(dimension, "bbob-largescale");
  nb_swaps = coco_get_nb_swaps(dimension, "bbob-largescale");*/

  coco_matrix_t B = coco_block_matrix(dimension, block_sizes.data(), nb_blocks);
  coco_compute_blockrotation(B.rows(), rseed, dimension, block_sizes.data(), nb_blocks);
  /*P1 = coco_allocate_vector_size_t(dimension);
  P2 = coco_allocate_vector_size_t(dimension);
  coco_compute_truncated_uniform_swap_permutation(P1, rseed + 2000000,
  dimension, nb_swaps, swap_range);
  coco_compute_truncated_uniform_swap_permutation(P2, rseed + 3000000,
  dimension, nb_swaps, swap_range);*/

  std::unique_ptr<f_gallagher_ls_t> gallagher = f_gallagher_problem_allocate(dimension, number_of_peaks);
  shared = &gallagher->shared; /* shortcut */
  coco_problem_ptr_t problem = std::move(gallagher);
  /* set shared fields needed later in transform_vars_gallagher_blockrotation.cpp */
  /* TODO: block-rotation related fields should be set in a separate function, not the definition of the problem*/
  idx_blocksize = 0;
  next_bs_change = block_sizes[idx_blocksize];
  for (i = 0; i < dimension; i++) {
    if (i >= next_bs_change) {
      idx_blocksize++;
      next_bs_change += block_sizes[idx_blocksize];
    }
    current_blocksize = block_sizes[idx_blocksize];
    shared->block_size_map[i] = current_blocksize;
    shared->first_non_zero_map[i] = next_bs_change - current_blocksize;
  }
  shared->B =
      coco_copy_block_matrix(B.rows(), dimension, block_sizes.data(), nb_blocks);

  rotation_problem = coco_problem_allocate(dimension, 1, 0);
  coco_problem_init_from_scalars(rotation_problem.get(), "dummy rotation", -5.0, 5.0, 0.0);
  rotation_problem =
      transform_vars_blockrotation(std::move(rotation_problem), B.rows(), dimension, block_sizes.data(), nb_blocks);

  std::vector<double> alpha_i_vals(number_of_peaks - 1);
  for (i = 0; i < number_of_peaks - 1; i++) {
    alpha_i_vals[i] = pow(1000, 2 * (double)(i) / ((double)(number_of_peaks - 2)));
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
  for (peak_index = 0; peak_index < number_of_peaks; peak_index++) {
    problem_i = &(shared->sub_problems[peak_index]);

    /* compute transformation parameters: */
    /* y_i and block-rotate it once and for all */
    std::vector<double> y_i(dimension);
    for (i = 0; i < dimension; i++) {
      y_i[i] = b * tmp_uniform[i + peak_index * dimension] - c;
      if (peak_index == 0) {
        y_i[i] *= a;
        problem->best_parameter[i] = 0;
        best_param_before_rotation[i] = y_i[i];
      }
    }
    transform_vars_blockrotation_apply(rotation_problem.get(), y_i.data(), y_i.data());
    if (peak_index == 0) {
      for (i = 0; i < dimension; i++) {
        best_param_after_rotation[i] = y_i[i];
      }
    }

    if (peak_index == 0) {
      alpha_i = first_condition;
    } else {
      alpha_i = alpha_i_vals[P_alpha_i[peak_index - 1]]; /*already square-rooted */
    }

    (*problem_i)->best_value[0] = 0; /* to prevent raising the assert */
    /* P_Lambda the permutation */
    std::vector<size_t> P_Lambda(dimension); /* random permutation of the values in C_i */
    bbob2009_unif(tmp_uniform2.data(), dimension, rseed + (long)(1000 * peak_index));
    coco_compute_permutation_from_sequence(P_Lambda.data(), tmp_uniform2.data(), dimension);
    /*coco_compute_random_permutation(P_Lambda, rseed + (long) (1000 *
     * peak_index), dimension);*/

    /* apply var transformations to sub problem*/
    *problem_i = transform_vars_scale(std::move(*problem_i), 1. / sqrt(sqrt(sqrt(alpha_i)))); /* sqrt( alpha^1/4) */
    *problem_i = transform_vars_conditioning(std::move(*problem_i), sqrt(alpha_i));
    /**problem_i = transform_vars_blockrotation(std::move(*problem_i), B_const, dimension,
     * block_sizes, nb_blocks);*/
    *problem_i = transform_vars_inverse_permutation(std::move(*problem_i), P_Lambda.data(), dimension);
    *problem_i = transform_vars_shift(std::move(*problem_i), y_i.data(), 0);
  }

  problem->evaluate_function(best_param_after_rotation.data(), problem->best_value.data());
  problem = transform_vars_blockrotation(std::move(problem), B.rows(), dimension, block_sizes.data(), nb_blocks);
  for (i = 0; i < dimension; i++) {
    problem->best_parameter[i] = best_param_before_rotation[i];
  }

  /*transform global objective function*/
  problem = transform_obj_oscillate(std::move(problem));
  problem = transform_obj_power(std::move(problem), 2.0);
  problem = transform_obj_penalize(std::move(problem), penalty_factor);
  problem = transform_obj_shift(std::move(problem), fopt);
  /*problem = transform_vars_gallagher_blockrotation(std::move(problem), shared);*/ /* block-matrix in shared data*/

  coco_problem_set_id(problem.get(), problem_id_template, function, instance, dimension);
  coco_problem_set_name(problem.get(), problem_name_template, function, instance, dimension);
  coco_problem_set_type(problem.get(), "5-weakly-structured");

  /*coco_free_memory(P1);
  coco_free_memory(P2);*/
  return problem;
}
