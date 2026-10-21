/**
 * @file suite_noisy_utilities.cpp
 * @brief Implementation of some functions (mostly handling instances) used by the bi-objective suites.
 * These are used throughout the COCO code base but should not be used by any external code.
 */

#include "coco.h"
#include "coco_internal.h"
#include "coco_random.cpp"
#include "coco_problem.cpp"
#include "suite_bbob_legacy_code.cpp"

namespace coco {

/***********************************************************************************************************/
/**
 * @name Methods and global variables needed to perform random sampling
 */
/**{@*/
static long RANDSEED = 30;  /** < @brief Random seed for sampling Uniform noise*/
static long RANDNSEED = 30; /** < @brief Random seed for sampling Gaussian noise*/

/**
 * @brief Increases the normal random seed by one unit
 * Needed to sample different values at each time
 */
void IncreaseRandomNSeed() {
  RANDNSEED += 1;
  if (RANDNSEED > static_cast<long>(1.0e9))
    RANDNSEED = 1;
}

/**
 * @brief Increases the uniform random seed by one unit
 * Needed to sample different values at each time
 */
void IncreaseRandomSeed() {
  RANDSEED += 1;
  if (RANDSEED > static_cast<long>(1.0e9))
    RANDSEED = 1;
}

/**
 * @brief Returns a sample from the gaussian distribution
 * using the legacy code for random number generations
 * Returns a double
 */
double SampleGaussianNoise() {
  double gaussian_noise_ptr[1] = {0.0};
  IncreaseRandomNSeed();
  Bbob2009Gauss(gaussian_noise_ptr, RANDNSEED);
  double gaussian_noise = gaussian_noise_ptr[0];
  return gaussian_noise;
}

/**
 * @brief Returns a sample from the uniform distribution
 * using the legacy code for random number generations
 * Returns a double
 */
double SampleUniformNoise() {
  double noise_vector[1] = {0.0};
  IncreaseRandomSeed();
  Bbob2009Unif(noise_vector, RANDSEED);
  double uniform_noise_term = noise_vector[0];
  return uniform_noise_term;
}
/**@}*/

/***********************************************************************************************************/
/**
 * @name Methods regarding boundary handling
 */

/**
 * @brief Applies a penalty to solution outside the feasible hypercube
 */
/**@{*/
double Problem::NoisyBoundaryPenalty(std::span<const double> x) const {
  double penalty, lower_bound, upper_bound;
  penalty = 0.0;
  for (size_t dimension = 0; dimension < number_of_variables_; dimension++) {
    lower_bound = smallest_values_of_interest_[dimension];
    upper_bound = largest_values_of_interest_[dimension];
    if (x[dimension] > upper_bound) {
      penalty += (x[dimension] - upper_bound) * (x[dimension] - upper_bound);
    } else if (x[dimension] < lower_bound) {
      penalty += (x[dimension] - lower_bound) * (x[dimension] - lower_bound);
    }
  }
  return 100.0 * penalty;
}
/**@}*/

/***********************************************************************************************************/

} // namespace coco

/* C interface (declared in coco.h) */

/**
 * @brief Resets both random seeds to the initial values
 */
void coco_reset_seeds() {
  coco::RANDSEED = 30;
  coco::RANDNSEED = 30;
}
