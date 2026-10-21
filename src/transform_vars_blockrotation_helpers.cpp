/**
 * @file transform_vars_blockrotation_helpers.cpp
 * @brief Implements the block rotations of the large-scale suite.
 */

#include <stdio.h>
#include <assert.h>
#include "coco.h"

#include "coco_random.cpp" /*tmp*/
#include "coco_utilities.cpp"
#include "suite_bbob_legacy_code.cpp" /*tmp*/

#include <time.h> /*tmp*/

namespace coco {

/* TODO: Document this file in doxygen style! */

/**
 * @brief Returns block size for block rotation matrices.
 *
 */
static size_t RotationMatrixBlockSize(size_t const dimension) {
  const double BLOCK_SIZE_RELATIVE = 1;      /* for block rotations, relative to dimension */
  const size_t MAX_BLOCK_SIZE_ABSOLUTE = 40; /* for block rotations */

  return DoubleToSizeT(
      DoubleMin(BLOCK_SIZE_RELATIVE * static_cast<double>(dimension), static_cast<double>(MAX_BLOCK_SIZE_ABSOLUTE)));
}

/**
 * @brief Computes a block-diagonal rotation matrix with blocks of sizes block_sizes, each a random
 * rotation based on seed.
 */
static BlockDiagonalMatrix ComputeBlockrotation(long seed, std::span<const size_t> block_sizes) {
  size_t current_blocksize, cumsum_prev_block_sizes;
  BlockDiagonalMatrix B(block_sizes);

  cumsum_prev_block_sizes = 0; /* shift in rows to account for the previous blocks */
  for (size_t idx_block = 0; idx_block < block_sizes.size(); idx_block++) {
    current_blocksize = block_sizes[idx_block];
    const Matrix current_block =
        Bbob2009ComputeRotation(seed + static_cast<long>(1000000) * static_cast<long>(idx_block), current_blocksize);

    /* now fill the block matrix*/
    for (size_t i = 0; i < current_blocksize; i++) {
      for (size_t j = 0; j < current_blocksize; j++) {
        B.InBlock(i + cumsum_prev_block_sizes, j) = current_block(i, j);
      }
    }

    cumsum_prev_block_sizes += current_blocksize;
  }
  return B;
}

/**
 * @brief returns the list of block_sizes and sets nb_blocks to its correct value
 *
 * Returns an empty list for unknown suites.
 */
static std::vector<size_t> GetBlockSizes(size_t *nb_blocks, size_t dimension, std::string_view suite_name) {
  std::vector<size_t> block_sizes;
  size_t block_size;

  if (suite_name == "bbob-largescale") {
    block_size = RotationMatrixBlockSize(dimension);
    *nb_blocks = dimension / block_size + ((dimension % block_size) > 0);
    block_sizes.resize(*nb_blocks);
    for (size_t i = 0; i < *nb_blocks - 1; i++) {
      block_sizes[i] = block_size;
    }
    block_sizes[*nb_blocks - 1] = dimension - (*nb_blocks - 1) * block_size; /*add rest*/
  }
  return block_sizes;
}

} // namespace coco
