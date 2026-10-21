/**
 * @file mo_utilities.cpp
 * @brief Definitions of miscellaneous functions used for multi-objective problems.
 */

#include <stdlib.h>
#include <stdio.h>
#include "coco.h"

/**
 * @brief Precision used when comparing multi-objective solutions.
 *
 * Two solutions are considered equal in objective space when their normalized difference is smaller than
 * mo_precision.
 *
 * @note mo_precision needs to be smaller than mo_discretization
 */
static const double mo_precision = 1e-13;

/**
 * @brief Discretization interval used for rounding normalized multi-objective solutions.
 *
 * @note mo_discretization needs to be larger than mo_precision
 */
static const double mo_discretization = 5 * 1e-13;

/**
 * @brief Computes and returns the Euclidean distance of the points first and second.
 */
static double mo_get_norm(std::span<const double> first, std::span<const double> second) {

  double norm = 0;

  assert(first.size() == second.size());
  for (size_t i = 0; i < first.size(); i++) {
    norm += pow(first[i] - second[i], 2);
  }

  return sqrt(norm);
}

/**
 * @brief Creates a rounded normalized version of the given solution w.r.t. the given ROI.
 *
 * If the solution seems to be better than the extremes it is corrected (2 objectives are assumed).
 */
static std::vector<double> mo_normalize(std::span<const double> y, std::span<const double> ideal,
                                        std::span<const double> nadir) {
  const size_t num_obj = y.size();
  assert(ideal.size() == num_obj && nadir.size() == num_obj);
  std::vector<double> normalized_y(num_obj);

  for (size_t i = 0; i < num_obj; i++) {
    assert((nadir[i] - ideal[i]) > mo_discretization);
    normalized_y[i] = (y[i] - ideal[i]) / (nadir[i] - ideal[i]);
    normalized_y[i] = coco_double_round(normalized_y[i] / mo_discretization) * mo_discretization;
    if (normalized_y[i] < 0) {
      coco_warning("mo_normalize(): Adjusting %.15e to %.15e", y[i], ideal[i]);
      normalized_y[i] = 0;
    }
  }

  for (size_t i = 0; i < num_obj; i++) {
    assert(num_obj == 2);
    if (coco_double_almost_equal(normalized_y[i], 0, mo_precision) && (normalized_y[1 - i] < 1)) {
      coco_warning("mo_normalize(): Adjusting %.15e to %.15e", y[1 - i], nadir[1 - i]);
      normalized_y[1 - i] = 1;
    }
  }

  return normalized_y;
}

/**
 * @brief Checks the dominance relation in the unconstrained minimization case between two normalized
 * solutions in the objective space.
 *
 * If two values are closer together than mo_precision, they are treated as equal.
 *
 * @return
 *  1 if normalized_y1 dominates normalized_y2 <br>
 *  0 if normalized_y1 and normalized_y2 are non-dominated <br>
 * -1 if normalized_y2 dominates normalized_y1 <br>
 * -2 if normalized_y1 is identical to normalized_y2
 */
static int mo_get_dominance(std::span<const double> normalized_y1, std::span<const double> normalized_y2) {

  int flag1 = 0;
  int flag2 = 0;

  assert(normalized_y1.size() == normalized_y2.size());
  for (size_t i = 0; i < normalized_y1.size(); i++) {
    if (coco_double_almost_equal(normalized_y1[i], normalized_y2[i], mo_precision)) {
      continue;
    } else if (normalized_y1[i] < normalized_y2[i]) {
      flag1 = 1;
    } else if (normalized_y1[i] > normalized_y2[i]) {
      flag2 = 1;
    }
  }

  if (flag1 && !flag2) {
    return 1;
  } else if (!flag1 && flag2) {
    return -1;
  } else if (flag1 && flag2) {
    return 0;
  } else { /* (!flag1 && !flag2) */
    return -2;
  }
}

/**
 * @brief Checks whether the normalized solution is within [0, 1]^num_obj.
 */
static int mo_is_within_ROI(std::span<const double> normalized_y) {

  int within = 1;

  for (size_t i = 0; i < normalized_y.size(); i++) {
    if (coco_double_almost_equal(normalized_y[i], 0, mo_precision) ||
        coco_double_almost_equal(normalized_y[i], 1, mo_precision) || (normalized_y[i] > 0 && normalized_y[i] < 1))
      continue;
    else
      within = 0;
  }
  return within;
}

/**
 * @brief Computes and returns the minimal normalized distance of the point normalized_y from the ROI
 * (equals 0 if within the ROI).
 *
 *  @note Assumes num_obj = 2 and normalized_y >= 0
 */
static double mo_get_distance_to_ROI(std::span<const double> normalized_y) {

  double diff_0, diff_1;

  if (mo_is_within_ROI(normalized_y))
    return 0;

  assert(normalized_y.size() == 2);
  assert(normalized_y[0] >= 0);
  assert(normalized_y[1] >= 0);

  diff_0 = normalized_y[0] - 1;
  diff_1 = normalized_y[1] - 1;
  if ((diff_0 > 0) && (diff_1 > 0)) {
    return sqrt(pow(diff_0, 2) + pow(diff_1, 2));
  } else if (diff_0 > 0)
    return diff_0;
  else
    return diff_1;
}
