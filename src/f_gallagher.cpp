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

namespace coco {

/**
 * @brief A random permutation type for the Gallagher problem.
 *
 * Needed to create a random permutation that is compatible with the old code.
 */
struct GallagherPermutation {
  double value;
  size_t index;
};

/**
 * Comparison function used for sorting.
 */
static int GallagherCompareDoubles(const void *a, const void *b) {
  double temp =
      static_cast<const GallagherPermutation *>(a)->value - static_cast<const GallagherPermutation *>(b)->value;
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
class BbobGallagher final : public BbobProblem<BbobGallagher> {
public:
  BbobGallagher(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                const BbobVariant variant, const size_t number_of_peaks, const double penalty_scale,
                const char *problem_id_template, const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "5-weakly-structured", problem_id_template, problem_name_template),
        number_of_peaks_(number_of_peaks), penalty_scale_(penalty_scale),
        fopt_(Bbob2009ComputeFopt(function, instance)) {
    /* Allocate temporary storage and space for the peaks */
    x_local_ = Matrix(number_of_peaks, dimension);
    arr_scales_ = Matrix(number_of_peaks, dimension);
    peak_values_.resize(number_of_peaks);
    tmx_.resize(dimension);

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
    std::vector<GallagherPermutation> rperm;

    if (number_of_peaks == peaks_101) {
      maxcondition1 = sqrt(maxcondition1);
      b = 10.;
      c = 5.;
    } else if (number_of_peaks == peaks_21) {
      b = 9.8;
      c = 4.9;
    } else {
      coco_error("BbobGallagher(): '%zu' is a non-supported number of peaks", number_of_peaks);
    }
    rotation_ = Bbob2009ComputeRotation(rseed, dimension);

    /* Initialize all the data of the inner problem */
    std::vector<double> random_numbers(number_of_peaks * dimension); /* This is large enough for all cases below */
    Bbob2009Unif(std::span(random_numbers).first(number_of_peaks - 1), rseed);
    rperm.resize(number_of_peaks - 1);
    for (size_t i = 0; i < number_of_peaks - 1; ++i) {
      rperm[i].value = random_numbers[i];
      rperm[i].index = i;
    }
    qsort(rperm.data(), number_of_peaks - 1, sizeof(rperm[0]), GallagherCompareDoubles);

    /* Random permutation */
    std::vector<double> arrCondition(number_of_peaks);
    arrCondition[0] = maxcondition1;
    peak_values_[0] = 10;
    for (size_t i = 1; i < number_of_peaks; ++i) {
      arrCondition[i] =
          pow(maxcondition, static_cast<double>(rperm[i - 1].index) / (static_cast<double>(number_of_peaks - 2)));
      peak_values_[i] =
          static_cast<double>(i - 1) / static_cast<double>(number_of_peaks - 2) * (fitvalues[1] - fitvalues[0]) +
          fitvalues[0];
    }

    rperm.resize(dimension);
    for (size_t i = 0; i < number_of_peaks; ++i) {
      Bbob2009Unif(random_numbers, rseed + static_cast<long>(1000 * i));
      for (size_t j = 0; j < dimension; ++j) {
        rperm[j].value = random_numbers[j];
        rperm[j].index = j;
      }
      qsort(rperm.data(), dimension, sizeof(rperm[0]), GallagherCompareDoubles);
      for (size_t j = 0; j < dimension; ++j) {
        arr_scales_(i, j) = pow(arrCondition[i], /* Lambda^alpha_i from the doc */
                                (static_cast<double>(rperm[j].index)) / (static_cast<double>(dimension - 1)) - 0.5);
      }
    }

    Bbob2009Unif(std::span(random_numbers).first(dimension * number_of_peaks), rseed);

    /* SBOX-COST wants optima in [-5, 5], the other suites in [-4, 4]. */
    double scale_factor = 0.8;
    if (variant == BbobVariant::BOXED) {
      scale_factor = 1.0;
    }

    for (size_t i = 0; i < dimension; ++i) {
      best_parameter_[i] = scale_factor * (b * random_numbers[i] - c);
      for (size_t j = 0; j < number_of_peaks; ++j) {
        x_local_(j, i) = 0.;
        for (size_t k = 0; k < dimension; ++k) {
          x_local_(j, i) += rotation_(i, k) * (b * random_numbers[j * dimension + k] - c);
        }
        if (j == 0) {
          x_local_(j, i) *= scale_factor;
        }
      }
    }

    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty =
        penalty_scale_ * BoundaryPenalty(x, smallest_values_of_interest_, largest_values_of_interest_);
    const double fac = -0.5 / static_cast<double>(number_of_variables_);
    double f = 0.;

    rotation_.Multiply(x, tmx_);
    for (size_t i = 0; i < number_of_peaks_; ++i) {
      const std::span<const double> peak = x_local_.Row(i);
      const std::span<const double> scales = arr_scales_.Row(i);
      double tmp2 = 0.;
      for (size_t j = 0; j < number_of_variables_; ++j) {
        const double tmp = (tmx_[j] - peak[j]);
        tmp2 += scales[j] * tmp * tmp;
      }
      tmp2 = peak_values_[i] * exp(fac * tmp2);
      f = DoubleMax(f, tmp2);
    }

    f = OscillateObjective(10. - f);
    const double y = f * f + penalty;
    return y + fopt_;
  }

private:
  const size_t number_of_peaks_;
  const double penalty_scale_;
  const double fopt_;
  Matrix rotation_;
  Matrix x_local_; /**< @brief The location of peak i in row i (number_of_peaks x dimension) */
  Matrix arr_scales_;
  std::vector<double> peak_values_;
  std::vector<double> tmx_; /**< @brief The rotated decision values */
};

/**
 * @brief The Gallagher functions f21 and f22 of the bbob-largescale suite with 101 and 21 peaks and a permuted block
 * rotation.
 *
 * f(x) = T_osz(10 - max_i w_i exp(-1 / (2 n) f_sphere(s_i Λ^(sqrt α_i) P_i (B x - y_i))))^2 + f_pen(x) + fopt with
 * the block rotation B, the peaks y_i, the random permutations P_i, the conditions α_i and s_i = α_i^(-1 / 8).
 */
class BbobGallagherPermblockdiag final : public BbobProblem<BbobGallagherPermblockdiag> {
public:
  BbobGallagherPermblockdiag(const size_t function, const size_t dimension, const size_t instance, const long rseed,
                             const size_t number_of_peaks, const char *problem_id_template,
                             const char *problem_name_template)
      : BbobProblem(function, dimension, instance, "5-weakly-structured", problem_id_template, problem_name_template),
        rotation_(Blockrotation(rseed, dimension)), fopt_(Bbob2009ComputeFopt(function, instance)), z_(dimension),
        z1_(dimension), z2_(dimension) {
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
      coco_error("BbobGallagherPermblockdiag(): '%zu' is a non-supported number of peaks", number_of_peaks);
    }

    std::vector<double> alpha_i_vals(number_of_peaks - 1);
    for (size_t i = 0; i < number_of_peaks - 1; i++) {
      alpha_i_vals[i] = pow(1000, 2 * static_cast<double>(i) / (static_cast<double>(number_of_peaks - 2)));
    }
    std::vector<double> tmp_uniform(dimension * number_of_peaks);
    Bbob2009Unif(std::span(tmp_uniform).first(number_of_peaks - 1), rseed);
    /* random permutation of the alpha_i's to allow sampling without replacement*/
    std::vector<size_t> P_alpha_i(number_of_peaks - 1);
    ComputePermutationFromSequence(P_alpha_i, std::span(tmp_uniform).first(number_of_peaks - 1));
    Bbob2009Unif(std::span(tmp_uniform).first(dimension * number_of_peaks), rseed);
    std::vector<double> tmp_uniform2(dimension);
    peaks_.reserve(number_of_peaks);
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
        best_parameter_ = y_i_unrotated;
      rotation_.Apply(y_i_unrotated, y_i);

      if (peak_index == 0) {
        alpha_i = first_condition;
      } else {
        alpha_i = alpha_i_vals[P_alpha_i[peak_index - 1]]; /*already square-rooted */
      }

      /* P_Lambda the permutation */
      std::vector<size_t> P_Lambda(dimension); /* random permutation of the values in C_i */
      Bbob2009Unif(tmp_uniform2, rseed + static_cast<long>(1000 * peak_index));
      ComputePermutationFromSequence(P_Lambda, tmp_uniform2);

      peaks_.push_back({VarsShift(y_i), VarsPermutation(P_Lambda, true), VarsConditioning(sqrt(alpha_i), dimension),
                        VarsScale(1. / sqrt(sqrt(sqrt(alpha_i))), dimension)}); /* sqrt( alpha^1/4) */
    }

    /* The best parameter is the first peak before the block rotation. Computing it with the inverse of the rotation
     * would add rounding errors (about 1e-11), which move the nadir points of the bi-objective mixed-integer suite
     * away from their reference values. */
    SetBestValueAtBestParameter();
  }

  double operator()(std::span<const double> x) {
    const double penalty = BoundaryPenalty(x, smallest_values_of_interest_, largest_values_of_interest_);
    double maxf = DBL_MAX;

    rotation_.Apply(x, z_);
    for (size_t i = 0; i < peaks_.size(); i++) {
      peaks_[i].shift.Apply(z_, z1_);
      peaks_[i].permutation.Apply(z1_, z2_);
      peaks_[i].conditioning.Apply(z2_, z2_);
      peaks_[i].scale.Apply(z2_, z2_);
      double y = SphereRaw(z2_);
      double w_i;
      if (i == 0) {
        w_i = 10;
      } else {
        w_i = 1.1 + 8.0 * ((static_cast<double>(i) + 1) - 2.0) / ((static_cast<double>(peaks_.size())) - 2.0);
      }
      y = w_i * exp(-1.0 / (2.0 * (static_cast<double>(number_of_variables_))) * y);
      if (i == 0 || maxf < y) {
        maxf = y;
      }
    }
    double y = pow(OscillateObjective(10.0 - maxf), exponent_);
    y += PENALTY_FACTOR * penalty;
    return y + fopt_;
  }

private:
  /** @brief The transformations of a peak, which are applied in this order. */
  struct Peak {
    VarsShift shift; /**< @brief Subtracts the location of the peak */
    VarsPermutation permutation;
    VarsConditioning conditioning;
    VarsScale scale;
  };

  static BlockDiagonalMatrix Blockrotation(const long rseed, const size_t dimension) {
    size_t nb_blocks;
    const std::vector<size_t> block_sizes = GetBlockSizes(&nb_blocks, dimension, "bbob-largescale");
    return ComputeBlockrotation(rseed, block_sizes);
  }

  static constexpr double PENALTY_FACTOR = 1.0;
  /**
   * @brief The exponent of the objective value. It is not a constant because compilers replace pow(y, 2.0) by
   * y * y, which can differ from pow() in the last bit.
   */
  const double exponent_ = 2.0;

  const VarsBlockrotation rotation_;
  std::vector<Peak> peaks_;
  const double fopt_;
  std::vector<double> z_, z1_, z2_; /**< @brief The transformed decision values */
};

} // namespace coco
