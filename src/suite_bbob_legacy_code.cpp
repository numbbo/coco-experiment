/**
 * @file suite_bbob_legacy_code.cpp
 * @brief Legacy code from BBOB2009 required to replicate the 2009 functions.
 *
 * All of this code should only be used by the suite_bbob2009 functions to provide compatibility to the
 * legacy code. New test beds should strive to use the new COCO facilities for random number generation etc.
 */

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "coco.h"

/** @brief Computes the minimum of the two values. */
static double bbob2009_fmin(double a, double b) { return (a < b) ? a : b; }

/** @brief Computes the maximum of the two values. */
static double bbob2009_fmax(double a, double b) { return (a > b) ? a : b; }

/** @brief Rounds the given value. */
static double bbob2009_round(double x) { return floor(x + 0.5); }

/**
 * @brief Generates N uniform random numbers using inseed as the seed and stores them in r.
 */
static void bbob2009_unif(double *r, size_t N, long inseed) {
  /* generates N uniform numbers with starting seed */
  long tmp;
  long rgrand[32];

  if (inseed < 0)
    inseed = -inseed;
  if (inseed < 1)
    inseed = 1;
  long aktseed = inseed;
  for (long i = 39; i >= 0; i--) {
    tmp = static_cast<int>(floor(static_cast<double>(aktseed) / static_cast<double>(127773)));
    aktseed = 16807 * (aktseed - tmp * 127773) - 2836 * tmp;
    if (aktseed < 0)
      aktseed = aktseed + 2147483647;
    if (i < 32)
      rgrand[i] = aktseed;
  }
  long aktrand = rgrand[0];
  for (long i = 0; i < static_cast<long>(N); i++) {
    tmp = static_cast<int>(floor(static_cast<double>(aktseed) / static_cast<double>(127773)));
    aktseed = 16807 * (aktseed - tmp * 127773) - 2836 * tmp;
    if (aktseed < 0)
      aktseed = aktseed + 2147483647;
    tmp = static_cast<int>(floor(static_cast<double>(aktrand) / static_cast<double>(67108865)));
    aktrand = rgrand[tmp];
    rgrand[tmp] = aktseed;
    r[i] = static_cast<double>(aktrand) / 2.147483647e9;
    if (r[i] == 0.) {
      r[i] = 1e-99;
    }
  }
  return;
}

/**
 * @brief Converts from packed (column-major) storage to a matrix.
 */
static void bbob2009_reshape(coco_matrix_t &B, const double *vector, const size_t m, const size_t n) {
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      B(i, j) = vector[j * m + i];
    }
  }
}

/**
 * @brief Generates N Gaussian random numbers using the given seed and stores them in g.
 */
static void bbob2009_gauss(double *g, const size_t N, const long seed) {
  std::vector<double> uniftmp(2 * N);
  bbob2009_unif(uniftmp.data(), 2 * N, seed);

  for (size_t i = 0; i < N; i++) {
    g[i] = sqrt(-2 * log(uniftmp[i])) * cos(2 * coco_pi * uniftmp[N + i]);
    if (g[i] == 0.)
      g[i] = 1e-99;
  }
  return;
}

/**
 * @brief Computes a DIM by DIM rotation matrix based on seed.
 */
static coco_matrix_t bbob2009_compute_rotation(const long seed, const size_t DIM) {
  coco_matrix_t B(DIM, DIM);
  double prod;
  std::vector<double> gvect(DIM * DIM);

  bbob2009_gauss(gvect.data(), DIM * DIM, seed);
  bbob2009_reshape(B, gvect.data(), DIM, DIM);
  /*1st coordinate is row, 2nd is column.*/

  for (size_t i = 0; i < DIM; i++) {
    for (size_t j = 0; j < i; j++) {
      prod = 0;
      for (size_t k = 0; k < DIM; k++)
        prod += B(k, i) * B(k, j);
      for (size_t k = 0; k < DIM; k++)
        B(k, i) -= prod * B(k, j);
    }
    prod = 0;
    for (size_t k = 0; k < DIM; k++)
      prod += B(k, i) * B(k, i);
    for (size_t k = 0; k < DIM; k++)
      B(k, i) /= sqrt(prod);
  }
  return B;
}

/**
 * @brief Randomly computes the location of the global optimum.
 */
static void bbob2009_compute_xopt(double *xopt, const long seed, const size_t DIM) {
  bbob2009_unif(xopt, DIM, seed);
  for (size_t i = 0; i < DIM; i++) {
    xopt[i] = 8 * floor(1e4 * xopt[i]) / 1e4 - 4;
    if (xopt[i] == 0.0)
      xopt[i] = -1e-5;
  }
}

/**
 * @brief Randomly computes the location of the global optimum in the larger space
 * of [-5, 5]^DIM compared to `bbob2009_compute_xopt`.
 */
static void sbox_cost_compute_xopt(double *xopt, const long seed, const size_t DIM) {
  bbob2009_unif(xopt, DIM, seed);
  for (size_t i = 0; i < DIM; i++) {
    xopt[i] = 10 * floor(1e4 * xopt[i]) / 1e4 - 5;
    if (xopt[i] == 0.0)
      xopt[i] = -1e-5;
  }
}

/**
 * @brief The suite a bbob function is created for.
 *
 * The functions differ slightly between the suites: bbob-boxed places the optimum in [-5, 5]^DIM and
 * bbob-largescale normalizes the objective value by the dimension and scales by the block size.
 */
enum class coco_bbob_variant_t { bbob, boxed, largescale };

/**
 * @brief Randomly computes the location of the global optimum for the given suite variant.
 */
static void bbob_compute_xopt(double *xopt, const long seed, const size_t DIM, const coco_bbob_variant_t variant) {
  if (variant == coco_bbob_variant_t::boxed)
    sbox_cost_compute_xopt(xopt, seed, DIM);
  else
    bbob2009_compute_xopt(xopt, seed, DIM);
}

/**
 * @brief Randomly chooses the objective offset for the given function and instance.
 */
static double bbob2009_compute_fopt(const size_t function, const size_t instance) {
  long rseed, rrseed;
  double gval, gval2;

  if (function == 4)
    rseed = 3;
  else if (function == 18)
    rseed = 17;
  else if (function == 101 || function == 102 || function == 103 || function == 107 || function == 108 ||
           function == 109)
    rseed = 1;
  else if (function == 104 || function == 105 || function == 106 || function == 110 || function == 111 ||
           function == 112)
    rseed = 8;
  else if (function == 113 || function == 114 || function == 115)
    rseed = 7;
  else if (function == 116 || function == 117 || function == 118)
    rseed = 10;
  else if (function == 119 || function == 120 || function == 121)
    rseed = 14;
  else if (function == 122 || function == 123 || function == 124)
    rseed = 17;
  else if (function == 125 || function == 126 || function == 127)
    rseed = 19;
  else if (function == 128 || function == 129 || function == 130)
    rseed = 21;
  else
    rseed = static_cast<long>(function);

  rrseed = rseed + static_cast<long>(10000 * instance);
  bbob2009_gauss(&gval, 1, rrseed);
  bbob2009_gauss(&gval2, 1, rrseed + 1);
  return bbob2009_fmin(1000., bbob2009_fmax(-1000., bbob2009_round(100. * 100. * gval / gval2) / 100.));
}
