#pragma once
#include "coco.h"
#undef INTERFACE
size_t *coco_get_block_sizes(size_t *nb_blocks, size_t dimension, const char *suite_name);
double **coco_copy_block_matrix(const double *const *B, const size_t dimension, const size_t *block_sizes,
                                const size_t nb_blocks);
void coco_compute_blockrotation(double **B, long seed, COCO_UNUSED size_t n, size_t *block_sizes, size_t nb_blocks);
void coco_free_block_matrix(double **matrix, const size_t n);
double **coco_allocate_blockmatrix(const size_t n, const size_t *block_sizes, const size_t nb_blocks);
size_t coco_rotation_matrix_block_size(size_t const dimension);
