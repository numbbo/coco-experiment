#pragma once
#include "coco.h"


size_t coco_get_nb_swaps(size_t dimension, const char *suite_name);
size_t coco_get_swap_range(size_t dimension, const char *suite_name);
size_t *coco_duplicate_size_t_vector(const size_t *src, const size_t number_of_elements);
void coco_compute_truncated_uniform_swap_permutation(size_t *P, long seed, size_t n, size_t nb_swaps,
                                                     size_t swap_range);
long coco_random_unif_integer(long lower_bound, long upper_bound, long seed);
void coco_compute_permutation_from_sequence(size_t *P, double *seq, size_t length);
void coco_compute_random_permutation(size_t *P, long seed, size_t n);
// FIXME: Remove unneeded forward declarations
// int f_compare_doubles_for_random_permutation(const void *a, const void *b);
