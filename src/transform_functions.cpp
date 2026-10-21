/**
 * @file transform_functions.cpp
 * @brief Transformations of the decision values as plain function objects.
 *
 * Each transformation computes y = T(x) with apply() and x = T^-1(y) with inverse(). The element-wise ones
 * allow x and y to be the same array. The bbob problem classes apply them explicitly; the transform_vars_*
 * problems wrap them around an inner problem.
 */

#include <math.h>

#include "coco.h"
#include "brentq.cpp"
#include "coco_utilities.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_vars_blockrotation_helpers.cpp"
#include "transform_vars_permutation_helpers.cpp"

/**
 * @brief The oscillating transformation T_osz of a single value.
 */
static double tosz_uv(const double xi, const double alpha) {
  double yi;
  double tmp, base;
  if (xi > 0.0) {
    tmp = log(xi) / alpha;
    base = exp(tmp + 0.49 * (sin(tmp) + sin(0.79 * tmp)));
    yi = pow(base, alpha);
  } else if (xi < 0.0) {
    tmp = log(-xi) / alpha;
    base = exp(tmp + 0.49 * (sin(0.55 * tmp) + sin(0.31 * tmp)));
    yi = -pow(base, alpha);
  } else {
    yi = 0.0;
  }
  return yi;
}

/**
 * @brief The asymmetric transformation T_asy of the value xi with index i of n values.
 */
static double tasy_uv(const double xi, const double beta, const size_t i, const size_t n) {
  double yi;
  double exponent;
  if (xi > 0.0) {
    exponent = 1.0 + ((beta * static_cast<double>(i)) / (static_cast<double>(n) - 1.0)) * sqrt(xi);
    yi = pow(xi, exponent);
  } else {
    yi = xi;
  }
  return yi;
}

/**
 * @brief Shifts the decision values: y = x - offset.
 */
class vars_shift_t {
public:
  explicit vars_shift_t(const double *offset, const size_t n) : offset(offset, offset + n) {}

  void apply(const double *x, double *y) const {
    for (size_t i = 0; i < offset.size(); ++i)
      y[i] = x[i] - offset[i];
  }

  void inverse(const double *y, double *x) const {
    for (size_t i = 0; i < offset.size(); ++i)
      x[i] = y[i] + offset[i];
  }

private:
  std::vector<double> offset;
};

/**
 * @brief Scales the decision values by alpha^(0.5 i / (n - 1)).
 */
class vars_conditioning_t {
public:
  vars_conditioning_t(const double alpha, const size_t n) : scales(n) {
    for (size_t i = 0; i < n; ++i)
      scales[i] = pow(alpha, 0.5 * static_cast<double>(i) / (static_cast<double>(n) - 1.0));
  }

  void apply(const double *x, double *y) const {
    for (size_t i = 0; i < scales.size(); ++i)
      y[i] = scales[i] * x[i];
  }

  void inverse(const double *y, double *x) const {
    for (size_t i = 0; i < scales.size(); ++i)
      x[i] = y[i] / scales[i];
  }

private:
  std::vector<double> scales; /**< @brief alpha^(0.5 i / (n - 1)) */
};

/**
 * @brief The oscillating transformation T_osz.
 */
class vars_oscillate_t {
public:
  explicit vars_oscillate_t(const size_t n) : n(n) {}

  void apply(const double *x, double *y) const {
    for (size_t i = 0; i < n; ++i)
      y[i] = tosz_uv(x[i], alpha);
  }

  /** @brief Inverts the transformation numerically. */
  void inverse(const double *y, double *x) const {
    for (size_t i = 0; i < n; ++i)
      x[i] = brentinv([](const double xi) { return tosz_uv(xi, alpha); }, y[i]);
  }

private:
  static constexpr double alpha = 0.1;
  size_t n;
};

/**
 * @brief The asymmetric transformation T_asy^beta.
 */
class vars_asymmetric_t {
public:
  vars_asymmetric_t(const double beta, const size_t n) : beta(beta), n(n) {}

  void apply(const double *x, double *y) const {
    for (size_t i = 0; i < n; ++i)
      y[i] = tasy_uv(x[i], beta, i, n);
  }

  /** @brief Inverts the transformation numerically. */
  void inverse(const double *y, double *x) const {
    for (size_t i = 0; i < n; ++i)
      x[i] = brentinv([this, i](const double xi) { return tasy_uv(xi, beta, i, n); }, y[i]);
  }

private:
  double beta;
  size_t n;
};

/**
 * @brief The affine transformation y = M x + b; x and y must not be the same array.
 */
class vars_affine_t {
public:
  vars_affine_t(coco_matrix_t M, std::vector<double> b) : M(std::move(M)), b(std::move(b)) {}

  void apply(const double *x, double *y) const { M.multiply_add(x, b.data(), y); }

  /** @brief Solves M x = y - b, which is M^T (y - b) if M is orthogonal. Returns false if M is singular. */
  bool inverse(const double *y, double *x) const {
    std::vector<double> shifted_y(b.size());
    for (size_t i = 0; i < b.size(); ++i)
      shifted_y[i] = y[i] - b[i];
    if (M.is_orthogonal()) {
      M.multiply_transposed(shifted_y.data(), x);
      return true;
    }
    return M.solve(shifted_y.data(), x);
  }

  const coco_matrix_t &matrix() const { return M; }

private:
  coco_matrix_t M; /**< @brief Has as many columns as x and as many rows as y */
  std::vector<double> b;
};

/**
 * @brief The block rotation y = B x; x and y must not be the same array.
 */
class vars_blockrotation_t {
public:
  explicit vars_blockrotation_t(coco_block_diagonal_matrix_t B) : B(std::move(B)) {}

  void apply(const double *x, double *y) const { B.multiply(x, y); }
  void inverse(const double *y, double *x) const { B.multiply_transposed(y, x); }

  const coco_block_diagonal_matrix_t &matrix() const { return B; }

private:
  coco_block_diagonal_matrix_t B;
};

/**
 * @brief Permutes the decision values: y[i] = x[P[i]], or y[P[i]] = x[i] for the inverse permutation; x and y
 * must not be the same array.
 */
class vars_permutation_t {
public:
  vars_permutation_t(const size_t *P, const size_t n, const bool is_inverse) : P(P, P + n), is_inverse(is_inverse) {}

  void apply(const double *x, double *y) const {
    if (is_inverse) {
      for (size_t i = 0; i < P.size(); ++i)
        y[P[i]] = x[i];
    } else {
      for (size_t i = 0; i < P.size(); ++i)
        y[i] = x[P[i]];
    }
  }

  void inverse(const double *y, double *x) const {
    for (size_t i = 0; i < P.size(); ++i) {
      if (is_inverse)
        x[i] = y[P[i]];
      else
        x[P[i]] = y[i];
    }
  }

private:
  std::vector<size_t> P; /**< @brief The permutation */
  bool is_inverse;
};

/**
 * @brief The permuted block rotation of the large-scale suite: y = P_second B P_first x, i.e. the decision values
 * are permuted, block rotated and permuted again; x and y must not be the same array.
 *
 * The permutations are truncated uniform swaps and the block sizes those of the bbob-largescale suite.
 */
class vars_permblockdiag_t {
public:
  vars_permblockdiag_t(const size_t n, const long seed_first, const long seed_rotation, const long seed_second)
      : first(random_permutation(n, seed_first).data(), n, false), rotation(random_blockrotation(n, seed_rotation)),
        second(random_permutation(n, seed_second).data(), n, false), tmp1(n), tmp2(n) {}

  void apply(const double *x, double *y) const {
    first.apply(x, tmp1.data());
    rotation.apply(tmp1.data(), tmp2.data());
    second.apply(tmp2.data(), y);
  }

  void inverse(const double *y, double *x) const {
    second.inverse(y, tmp2.data());
    rotation.inverse(tmp2.data(), tmp1.data());
    first.inverse(tmp1.data(), x);
  }

private:
  static std::vector<size_t> random_permutation(const size_t n, const long seed) {
    std::vector<size_t> P(n);
    coco_compute_truncated_uniform_swap_permutation(P.data(), seed, n, coco_get_nb_swaps(n, "bbob-largescale"),
                                                    coco_get_swap_range(n, "bbob-largescale"));
    return P;
  }

  static coco_block_diagonal_matrix_t random_blockrotation(const size_t n, const long seed) {
    size_t number_of_blocks;
    const std::vector<size_t> block_sizes = coco_get_block_sizes(&number_of_blocks, n, "bbob-largescale");
    return coco_compute_blockrotation(seed, block_sizes.data(), number_of_blocks);
  }

  vars_permutation_t first;
  vars_blockrotation_t rotation;
  vars_permutation_t second;
  mutable std::vector<double> tmp1, tmp2; /**< @brief Intermediate values */
};

/**
 * @brief Rounds the decision values for the step ellipsoid: to integers where |x| > 0.5 and to multiples of
 * 1 / alpha otherwise. Not invertible.
 */
class vars_round_step_t {
public:
  vars_round_step_t(const double alpha, const size_t n) : alpha(alpha), n(n) {}

  void apply(const double *x, double *y) const {
    for (size_t i = 0; i < n; ++i) {
      if (fabs(x[i]) > 0.5)
        y[i] = coco_double_round(x[i]);
      else
        y[i] = coco_double_round(alpha * x[i]) / alpha;
    }
  }

private:
  double alpha;
  size_t n;
};

/**
 * @brief The penalty of the large-scale suite for x outside the region of interest [lower, upper].
 */
static double coco_boundary_penalty(const double *x, const double *lower, const double *upper, const size_t n) {
  double penalty = 0.0;
  for (size_t i = 0; i < n; ++i) {
    const double c1 = x[i] - upper[i];
    const double c2 = lower[i] - x[i];
    assert(lower[i] < upper[i]);
    if (c1 > 0.0) {
      penalty += c1 * c1;
    } else if (c2 > 0.0) {
      penalty += c2 * c2;
    }
  }
  return penalty;
}

/**
 * @brief The factor with which the large-scale suite normalizes the objective values by the dimension.
 */
static double coco_norm_by_dim_factor(const size_t n) { return bbob2009_fmin(1, 40. / static_cast<double>(n)); }
