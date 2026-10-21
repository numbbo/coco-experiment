/**
 * @file f_gallagher.cpp
 * @brief Implementation of the Gallagher function and problem.
 */

#include <assert.h>
#include <math.h>

#include "coco.h"
#include "bbob_problem.cpp"
#include "f_sphere.cpp"
#include "transform_vars_blockrotation_helpers.cpp"
#include "transform_vars_permutation_helpers.cpp"

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
 * @brief The BBOB Gallagher functions f21 and f22 with 101 and 21 peaks, also in the bbob-boxed and bbob-noisy
 * suites.
 *
 * f(x) = T_osz(10 - max_i w_i exp(-1 / (2 n) (R x - y_i)^T C_i (R x - y_i)))^2 + penalty_scale f_pen(x) + fopt with
 * the peaks y_i, which lie in [-5, 5] for the boxed variant and in [-4, 4] otherwise.
 */
class bbob_gallagher_t final : public bbob_problem_t<bbob_gallagher_t> {
public:
  bbob_gallagher_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                   const coco_bbob_variant_t variant, const size_t number_of_peaks, const double penalty_scale,
                   const char *problem_id_template, const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "5-weakly-structured", problem_id_template,
                       problem_name_template),
        number_of_peaks(number_of_peaks), penalty_scale(penalty_scale),
        fopt(bbob2009_compute_fopt(function, instance)) {
    /* Allocate temporary storage and space for the peaks */
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
      coco_error("bbob_gallagher_t(): '%zu' is a non-supported number of peaks", number_of_peaks);
    }
    rotation = bbob2009_compute_rotation(rseed, dimension);

    /* Initialize all the data of the inner problem */
    std::vector<double> random_numbers(number_of_peaks * dimension); /* This is large enough for all cases below */
    bbob2009_unif(std::span(random_numbers).first(number_of_peaks - 1), rseed);
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
      bbob2009_unif(random_numbers, rseed + static_cast<long>(1000 * i));
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

    bbob2009_unif(std::span(random_numbers).first(dimension * number_of_peaks), rseed);

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

    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty =
        penalty_scale * coco_boundary_penalty(x, smallest_values_of_interest, largest_values_of_interest);
    const double fac = -0.5 / static_cast<double>(number_of_variables);
    double f = 0.;

    rotation.multiply(x, tmx);
    for (size_t i = 0; i < number_of_peaks; ++i) {
      const std::span<const double> peak = x_local.row(i);
      const std::span<const double> scales = arr_scales.row(i);
      double tmp2 = 0.;
      for (size_t j = 0; j < number_of_variables; ++j) {
        const double tmp = (tmx[j] - peak[j]);
        tmp2 += scales[j] * tmp * tmp;
      }
      tmp2 = peak_values[i] * exp(fac * tmp2);
      f = coco_double_max(f, tmp2);
    }

    f = coco_oscillate_objective(10. - f);
    const double y = f * f + penalty;
    return y + fopt;
  }

private:
  const size_t number_of_peaks;
  const double penalty_scale;
  const double fopt;
  coco_matrix_t rotation;
  coco_matrix_t x_local; /**< @brief The location of peak i in row i (number_of_peaks x dimension) */
  coco_matrix_t arr_scales;
  std::vector<double> peak_values;
  std::vector<double> tmx; /**< @brief The rotated decision values */
};

/**
 * @brief The Gallagher functions f21 and f22 of the bbob-largescale suite with 101 and 21 peaks and a permuted block
 * rotation.
 *
 * f(x) = T_osz(10 - max_i w_i exp(-1 / (2 n) f_sphere(s_i Λ^(sqrt α_i) P_i (B x - y_i))))^2 + f_pen(x) + fopt with
 * the block rotation B, the peaks y_i, the random permutations P_i, the conditions α_i and s_i = α_i^(-1 / 8).
 */
class bbob_gallagher_permblockdiag_t final : public bbob_problem_t<bbob_gallagher_permblockdiag_t> {
public:
  bbob_gallagher_permblockdiag_t(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                                 const size_t number_of_peaks, const char *problem_id_template,
                                 const char *problem_name_template)
      : bbob_problem_t(function, dimension, instance, "5-weakly-structured", problem_id_template,
                       problem_name_template),
        rotation(blockrotation(rseed, dimension)), fopt(bbob2009_compute_fopt(function, instance)), z(dimension),
        z1(dimension), z2(dimension) {
    const size_t peaks_21 = 21;
    const size_t peaks_101 = 101;
    double a = 0.8, b = 0, c = 0;
    double first_condition = 0;
    double alpha_i;
    /* first_condition satisfies the old code and the doc but seems wrong in that
     * it is, with very high probability, not the largest condition level!!! */

    if (number_of_peaks == peaks_101) {
      first_condition = 1000;
      b = 10.0;
      c = 5.0;
    } else if (number_of_peaks == peaks_21) {
      first_condition = 1000 * 1000;
      b = 9.8;
      c = 4.9;
    } else {
      coco_error("bbob_gallagher_permblockdiag_t(): '%zu' is a non-supported number of peaks", number_of_peaks);
    }

    std::vector<double> alpha_i_vals(number_of_peaks - 1);
    for (size_t i = 0; i < number_of_peaks - 1; i++) {
      alpha_i_vals[i] = pow(1000, 2 * static_cast<double>(i) / (static_cast<double>(number_of_peaks - 2)));
    }
    std::vector<double> tmp_uniform(dimension * number_of_peaks);
    bbob2009_unif(std::span(tmp_uniform).first(number_of_peaks - 1), rseed);
    /* random permutation of the alpha_i's to allow sampling without replacement*/
    std::vector<size_t> P_alpha_i(number_of_peaks - 1);
    coco_compute_permutation_from_sequence(P_alpha_i, std::span(tmp_uniform).first(number_of_peaks - 1));
    bbob2009_unif(std::span(tmp_uniform).first(dimension * number_of_peaks), rseed);
    std::vector<double> tmp_uniform2(dimension);
    peaks.reserve(number_of_peaks);
    for (size_t peak_index = 0; peak_index < number_of_peaks; peak_index++) {
      /* compute transformation parameters: */
      /* y_i and block-rotate it once and for all */
      std::vector<double> y_i(dimension);
      for (size_t i = 0; i < dimension; i++) {
        y_i[i] = b * tmp_uniform[i + peak_index * dimension] - c;
        if (peak_index == 0)
          y_i[i] *= a;
      }
      const std::vector<double> y_i_unrotated = y_i;
      if (peak_index == 0)
        best_parameter = y_i_unrotated;
      rotation.apply(y_i_unrotated, y_i);

      if (peak_index == 0) {
        alpha_i = first_condition;
      } else {
        alpha_i = alpha_i_vals[P_alpha_i[peak_index - 1]]; /*already square-rooted */
      }

      /* P_Lambda the permutation */
      std::vector<size_t> P_Lambda(dimension); /* random permutation of the values in C_i */
      bbob2009_unif(tmp_uniform2, rseed + static_cast<long>(1000 * peak_index));
      coco_compute_permutation_from_sequence(P_Lambda, tmp_uniform2);

      peaks.push_back({vars_shift_t(y_i), vars_permutation_t(P_Lambda, true),
                       vars_conditioning_t(sqrt(alpha_i), dimension),
                       vars_scale_t(1. / sqrt(sqrt(sqrt(alpha_i))), dimension)}); /* sqrt( alpha^1/4) */
    }

    /* The best parameter is the first peak before the block rotation. Computing it with the inverse of the rotation
     * would add rounding errors (about 1e-11), which move the nadir points of the bi-objective mixed-integer suite
     * away from their reference values. */
    set_best_value_at_best_parameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = coco_boundary_penalty(x, smallest_values_of_interest, largest_values_of_interest);
    double maxf = DBL_MAX;

    rotation.apply(x, z);
    for (size_t i = 0; i < peaks.size(); i++) {
      peaks[i].shift.apply(z, z1);
      peaks[i].permutation.apply(z1, z2);
      peaks[i].conditioning.apply(z2, z2);
      peaks[i].scale.apply(z2, z2);
      double y = f_sphere_raw(z2);
      double w_i;
      if (i == 0) {
        w_i = 10;
      } else {
        w_i = 1.1 + 8.0 * ((static_cast<double>(i) + 1) - 2.0) / ((static_cast<double>(peaks.size())) - 2.0);
      }
      y = w_i * exp(-1.0 / (2.0 * (static_cast<double>(number_of_variables))) * y);
      if (i == 0 || maxf < y) {
        maxf = y;
      }
    }
    double y = pow(coco_oscillate_objective(10.0 - maxf), exponent);
    y += penalty_factor * penalty;
    return y + fopt;
  }

private:
  /** @brief The transformations of a peak, which are applied in this order. */
  struct peak_t {
    vars_shift_t shift; /**< @brief Subtracts the location of the peak */
    vars_permutation_t permutation;
    vars_conditioning_t conditioning;
    vars_scale_t scale;
  };

  static coco_block_diagonal_matrix_t blockrotation(const long rseed, const size_t dimension) {
    size_t nb_blocks;
    const std::vector<size_t> block_sizes = coco_get_block_sizes(&nb_blocks, dimension, "bbob-largescale");
    return coco_compute_blockrotation(rseed, block_sizes);
  }

  static constexpr double penalty_factor = 1.0;
  /**
   * @brief The exponent of the objective value. It is not a constant because compilers replace pow(y, 2.0) by
   * y * y, which can differ from pow() in the last bit.
   */
  const double exponent = 2.0;

  const vars_blockrotation_t rotation;
  std::vector<peak_t> peaks;
  const double fopt;
  std::vector<double> z, z1, z2; /**< @brief The transformed decision values */
};
