/**
 * @file transform_vars_permutation_helpers.cpp
 * @brief Implements the random permutations of the large-scale suite.
 */

#include <stdio.h>
#include <assert.h>
#include "coco.h"

#include "coco_random.cpp"
#include "coco_utilities.cpp"
#include "suite_bbob_legacy_code.cpp" /*tmp*/

#include <time.h> /*tmp*/

#include <algorithm>

namespace coco {

/* TODO: Document this file in doxygen style! */

/**
 * @brief Sorts the indices in P by the corresponding values. Used as a random permutation generator when the
 * values are random. Indices with equal values keep their order.
 */
static void SortIndicesByValues(std::span<size_t> P, std::span<const double> values) {
  std::stable_sort(P.begin(), P.end(), [values](const size_t a, const size_t b) { return values[a] < values[b]; });
}

/**
 * @brief generates a random, uniformly sampled, permutation and puts it in P
 * Wassim: move to coco_utilities?
 */
static void ComputeRandomPermutation(std::span<size_t> P, long seed) {

  std::vector<double> random_data(P.size());
  Bbob2009Gauss(random_data, seed);
  for (size_t i = 0; i < P.size(); i++) {
    P[i] = i;
  }
  SortIndicesByValues(P, random_data);
}

/**
 * @brief generates a permutation by sorting a sequence and puts it in P
 */
static void ComputePermutationFromSequence(std::span<size_t> P, std::span<const double> seq) {
  assert(P.size() == seq.size());
  for (size_t i = 0; i < P.size(); i++) {
    P[i] = i;
  }
  SortIndicesByValues(P, seq);
}

/**
 * @brief returns a uniformly distributed integer between lower_bound and upper_bound using seed
 * without using coco_random_new.
 * Move to coco_utilities?
 */
static long RandomUnifInteger(long lower_bound, long upper_bound, long seed) {
  long range, rand_int;
  std::vector<double> tmp_uniform(1);
  Bbob2009Unif(tmp_uniform, seed);
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
static void ComputeTruncatedUniformSwapPermutation(std::span<size_t> P, long seed, size_t nb_swaps, size_t swap_range) {
  const size_t n = P.size();
  size_t lower_bound, upper_bound, first_swap_var, second_swap_var, tmp;

  if (n <= 40) {
    /* Do an identity permutation for dimensions <= 40 */
    for (size_t i = 0; i < n; i++)
      P[i] = i;
    return;
  }

  std::vector<double> random_data(n);
  Bbob2009Unif(random_data, seed);

  std::vector<size_t> idx_order(n);
  for (size_t i = 0; i < n; i++) {
    P[i] = i;
    idx_order[i] = i;
  }

  if (swap_range > 0) {
    /*sort the random data in random_data and arrange idx_order accordingly*/
    /*did not use ComputeRandomPermutation to only use the seed once*/
    SortIndicesByValues(idx_order, random_data);
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

      second_swap_var =
          static_cast<size_t>(RandomUnifInteger(static_cast<long>(lower_bound), static_cast<long>(upper_bound) - 1,
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
    ComputeRandomPermutation(P, seed);
  }
}

/**
 * @brief return the swap_range corresponding to the problem in the given suite
 */
static size_t GetSwapRange(size_t dimension, std::string_view suite_name) {
  if (suite_name == "bbob-largescale") {
    return dimension / 3;
  } else {
    coco_error("GetSwapRange(): unknown problem suite");
    return 0;
  }
}

/**
 * @brief return the number of swaps corresponding to the problem in the given suite
 */
size_t GetNbSwaps(size_t dimension, std::string_view suite_name) {
  if (suite_name == "bbob-largescale") {
    return dimension;
  } else {
    coco_error("GetNbSwaps(): unknown problem suite");
    return 0;
  }
}

} // namespace coco
