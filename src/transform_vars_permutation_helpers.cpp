/**
 * @file transform_vars_permutation_helpers.cpp
 * @brief implements functions needed by transform_vars_permutation.cpp
 */

#include <stdio.h>
#include <assert.h>
#include "coco.h"

#include "coco_random.cpp"
#include "coco_utilities.cpp"
#include "suite_bbob_legacy_code.cpp" /*tmp*/

#include <time.h> /*tmp*/

#include <algorithm>

/* TODO: Document this file in doxygen style! */

/**
 * @brief Sorts the indices in P (of length n) by the corresponding values. Used as a random permutation
 * generator when the values are random. Indices with equal values keep their order.
 */
static void coco_sort_indices_by_values(size_t *P, const double *values, const size_t n) {
  std::stable_sort(P, P + n, [values](const size_t a, const size_t b) { return values[a] < values[b]; });
}

/**
 * @brief generates a random, uniformly sampled, permutation and puts it in P
 * Wassim: move to coco_utilities?
 */
static void coco_compute_random_permutation(size_t *P, long seed, size_t n) {

  std::vector<double> random_data(n);
  bbob2009_gauss(random_data.data(), n, seed);
  for (size_t i = 0; i < n; i++) {
    P[i] = i;
  }
  coco_sort_indices_by_values(P, random_data.data(), n);
}

/**
 * @brief generates a permutation by sorting a sequence and puts it in P
 */
static void coco_compute_permutation_from_sequence(size_t *P, const double *seq, size_t length) {

  for (size_t i = 0; i < length; i++) {
    P[i] = i;
  }
  coco_sort_indices_by_values(P, seq, length);
}

/**
 * @brief returns a uniformly distributed integer between lower_bound and upper_bound using seed
 * without using coco_random_new.
 * Move to coco_utilities?
 */
static long coco_random_unif_integer(long lower_bound, long upper_bound, long seed) {
  long range, rand_int;
  std::vector<double> tmp_uniform(1);
  bbob2009_unif(tmp_uniform.data(), 1, seed);
  range = upper_bound - lower_bound + 1;
  rand_int = (static_cast<long>(tmp_uniform[0] * static_cast<double>(range))) + lower_bound;
  return rand_int;
}

/**
 * @brief generates a random permutation resulting from nb_swaps truncated uniform swaps of range swap_range
 * missing parameters: dynamic_not_static pool, seems empirically irrelevant
 * for now so dynamic is implemented (simple since no need for tracking indices
 * if swap_range is 0, a random uniform permutation is generated
 */
static void coco_compute_truncated_uniform_swap_permutation(size_t *P, long seed, size_t n, size_t nb_swaps,
                                                            size_t swap_range) {
  size_t lower_bound, upper_bound, first_swap_var, second_swap_var, tmp;

  if (n <= 40) {
    /* Do an identity permutation for dimensions <= 40 */
    for (size_t i = 0; i < n; i++)
      P[i] = i;
    return;
  }

  std::vector<double> random_data(n);
  bbob2009_unif(random_data.data(), n, seed);

  std::vector<size_t> idx_order(n);
  for (size_t i = 0; i < n; i++) {
    P[i] = i;
    idx_order[i] = i;
  }

  if (swap_range > 0) {
    /*sort the random data in random_data and arrange idx_order accordingly*/
    /*did not use coco_compute_random_permutation to only use the seed once*/
    coco_sort_indices_by_values(idx_order.data(), random_data.data(), n);
    for (size_t idx_swap = 0; idx_swap < nb_swaps; idx_swap++) {
      first_swap_var = idx_order[idx_swap];
      if (first_swap_var < swap_range) {
        lower_bound = 0;
      } else {
        lower_bound = first_swap_var - swap_range;
      }
      if (first_swap_var + swap_range > n - 1) {
        upper_bound = n - 1;
      } else {
        upper_bound = first_swap_var + swap_range;
      }

      second_swap_var = static_cast<size_t>(coco_random_unif_integer(static_cast<long>(lower_bound),
                                                                     static_cast<long>(upper_bound) - 1,
                                                                     seed + static_cast<long>(1 + idx_swap) * 1000));
      if (second_swap_var >= first_swap_var) {
        second_swap_var += 1;
      }
      /* swap*/
      tmp = P[first_swap_var];
      P[first_swap_var] = P[second_swap_var];
      P[second_swap_var] = tmp;
    }
  } else {
    /* generate random permutation instead */
    coco_compute_random_permutation(P, seed, n);
  }
}

/**
 * @brief return the swap_range corresponding to the problem in the given suite
 */
static size_t coco_get_swap_range(size_t dimension, std::string_view suite_name) {
  if (suite_name == "bbob-largescale") {
    return dimension / 3;
  } else {
    coco_error("coco_get_swap_range(): unknown problem suite");
    return 0;
  }
}

/**
 * @brief return the number of swaps corresponding to the problem in the given suite
 */
size_t coco_get_nb_swaps(size_t dimension, std::string_view suite_name) {
  if (suite_name == "bbob-largescale") {
    return dimension;
  } else {
    coco_error("coco_get_nb_swaps(): unknown problem suite");
    return 0;
  }
}
