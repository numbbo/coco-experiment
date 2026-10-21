/**
 * @file coco_random.cpp
 * @brief Definitions of functions regarding COCO random numbers.
 *
 * @note This file contains non-C89-standard types (such as uint32_t and uint64_t), which should
 * eventually be fixed.
 */

#include "coco_fp_contract.h"

#include <math.h>

#include "coco.h"
#include <stdio.h>

/** @brief The type declared in coco.h, an empty base class of RandomState (see coco_internal.h). */
struct coco_random_state_s {};

namespace coco {

#define COCO_NORMAL_POLAR /* Use polar transformation method */

static constexpr size_t COCO_SHORT_LAG = 273;
static constexpr size_t COCO_LONG_LAG = 607;

/**
 * @brief A structure containing the state of the COCO random generator.
 */
struct RandomState : public coco_random_state_s {
public:
  explicit RandomState(uint32_t seed);

  /** @brief Returns the next uniformly distributed random number in [0, 1). */
  double Uniform();

  /** @brief Returns an approximately N(0, 1) distributed random number. */
  double Normal();

private:
  void Generate();

  double x_[COCO_LONG_LAG];
  size_t index_;
};

/**
 * @brief A lagged Fibonacci random number generator.
 *
 * This generator is nice because it is reasonably small and directly generates double values. The chosen
 * lags (607 and 273) lead to a generator with a period in excess of 2^607-1.
 */
void RandomState::Generate() {
  for (size_t i = 0; i < COCO_SHORT_LAG; ++i) {
    double t = x_[i] + x_[i + (COCO_LONG_LAG - COCO_SHORT_LAG)];
    if (t >= 1.0)
      t -= 1.0;
    x_[i] = t;
  }
  for (size_t i = COCO_SHORT_LAG; i < COCO_LONG_LAG; ++i) {
    double t = x_[i] + x_[i - COCO_SHORT_LAG];
    if (t >= 1.0)
      t -= 1.0;
    x_[i] = t;
  }
  index_ = 0;
}

RandomState::RandomState(uint32_t seed) {
  /* Expand seed to fill initial state array. */
  for (size_t i = 0; i < COCO_LONG_LAG; ++i) {
    /* Uses uint64_t to silence the compiler ("shift count negative or too big, undefined behavior" warning) */
    x_[i] = (static_cast<double>(seed)) / static_cast<double>((static_cast<uint64_t>(1UL) << 32) - 1);
    /* Advance seed based on simple RNG from TAOCP */
    seed = static_cast<uint32_t>(1812433253UL) * (seed ^ (seed >> 30)) + (static_cast<uint32_t>(i) + 1);
  }
  index_ = 12;
}

double RandomState::Uniform() {
  /* If we have consumed all random numbers in our archive, it is time to run the actual generator for one
   * iteration to refill the state with 'LONG_LAG' new values. */
  if (index_ >= COCO_LONG_LAG)
    Generate();
  return x_[index_++];
}

/**
 * Instead of using the (expensive) polar method, we may cheat and abuse the central limit theorem. The sum
 * of 12 uniform random values has mean 6, variance 1 and is approximately normal. Subtract 6 and you get
 * an approximately N(0, 1) random number.
 */
double RandomState::Normal() {
  double result;
#ifdef COCO_NORMAL_POLAR
  const double u1 = Uniform();
  const double u2 = Uniform();
  result = sqrt(-2 * log(u1)) * cos(2 * coco_pi * u2);
#else
  result = 0.0;
  for (int i = 0; i < 12; ++i) {
    result += uniform();
  }
  result -= 6.0;
#endif
  return result;
}

} // namespace coco

/* C interface (declared in coco.h) */

/** @brief Converts a pointer of the C interface to the internal class. */
inline coco::RandomState *FromC(coco_random_state_t *state) { return static_cast<coco::RandomState *>(state); }

coco_random_state_t *coco_random_new(uint32_t seed) { return new coco::RandomState(seed); }

void coco_random_free(coco_random_state_t *state) { delete FromC(state); }

double coco_random_uniform(coco_random_state_t *state) { return FromC(state)->Uniform(); }

double coco_random_normal(coco_random_state_t *state) { return FromC(state)->Normal(); }
