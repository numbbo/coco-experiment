/**
 * @file transform_functions.cpp
 * @brief Transformations of the decision values as plain function objects.
 *
 * Each transformation computes y = T(x) with apply() and x = T^-1(y) with inverse(). x and y have as many
 * elements as the transformation has variables. The element-wise ones allow x and y to be the same array. The bbob
 * problem classes apply them explicitly; the transform_vars_* problems wrap them around an inner problem.
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
  explicit vars_shift_t(std::span<const double> offset) : offset(offset.begin(), offset.end()) {}

  void apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == offset.size() && y.size() == offset.size());
    for (size_t i = 0; i < offset.size(); ++i)
      y[i] = x[i] - offset[i];
  }

  void inverse(std::span<const double> y, std::span<double> x) const {
    assert(y.size() == offset.size() && x.size() == offset.size());
    for (size_t i = 0; i < offset.size(); ++i)
      x[i] = y[i] + offset[i];
  }

private:
  std::vector<double> offset;
};

/**
 * @brief Scales the decision values: y = factor x.
 */
class vars_scale_t {
public:
  vars_scale_t(const double factor, const size_t n) : factor(factor), n(n) {}

  void apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == n && y.size() == n);
    for (size_t i = 0; i < n; ++i)
      y[i] = factor * x[i];
  }

  /** @brief Divides by the factor. Returns false if the factor is 0. */
  bool inverse(std::span<const double> y, std::span<double> x) const {
    assert(y.size() == n && x.size() == n);
    if (factor == 0.)
      return false;
    for (size_t i = 0; i < n; ++i)
      x[i] = y[i] / factor;
    return true;
  }

private:
  double factor;
  size_t n;
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

  void apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == scales.size() && y.size() == scales.size());
    for (size_t i = 0; i < scales.size(); ++i)
      y[i] = scales[i] * x[i];
  }

  void inverse(std::span<const double> y, std::span<double> x) const {
    assert(y.size() == scales.size() && x.size() == scales.size());
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
  void apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == y.size());
    for (size_t i = 0; i < x.size(); ++i)
      y[i] = tosz_uv(x[i], alpha);
  }

  /** @brief Inverts the transformation numerically. */
  void inverse(std::span<const double> y, std::span<double> x) const {
    assert(y.size() == x.size());
    for (size_t i = 0; i < y.size(); ++i)
      x[i] = brentinv([](const double xi) { return tosz_uv(xi, alpha); }, y[i]);
  }

private:
  static constexpr double alpha = 0.1;
};

/**
 * @brief The asymmetric transformation T_asy^beta.
 */
class vars_asymmetric_t {
public:
  vars_asymmetric_t(const double beta, const size_t n) : beta(beta), n(n) {}

  void apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == n && y.size() == n);
    for (size_t i = 0; i < n; ++i)
      y[i] = tasy_uv(x[i], beta, i, n);
  }

  /** @brief Inverts the transformation numerically. */
  void inverse(std::span<const double> y, std::span<double> x) const {
    assert(y.size() == n && x.size() == n);
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

  void apply(std::span<const double> x, std::span<double> y) const { M.multiply_add(x, b, y); }

  /** @brief Solves M x = y - b, which is M^T (y - b) if M is orthogonal. Returns false if M is singular. */
  bool inverse(std::span<const double> y, std::span<double> x) const {
    assert(y.size() == b.size());
    std::vector<double> shifted_y(b.size());
    for (size_t i = 0; i < b.size(); ++i)
      shifted_y[i] = y[i] - b[i];
    if (M.is_orthogonal()) {
      M.multiply_transposed(shifted_y, x);
      return true;
    }
    return M.solve(shifted_y, x);
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

  void apply(std::span<const double> x, std::span<double> y) const { B.multiply(x, y); }
  void inverse(std::span<const double> y, std::span<double> x) const { B.multiply_transposed(y, x); }

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
  vars_permutation_t(std::span<const size_t> P, const bool is_inverse)
      : P(P.begin(), P.end()), is_inverse(is_inverse) {}

  void apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == P.size() && y.size() == P.size());
    if (is_inverse) {
      for (size_t i = 0; i < P.size(); ++i)
        y[P[i]] = x[i];
    } else {
      for (size_t i = 0; i < P.size(); ++i)
        y[i] = x[P[i]];
    }
  }

  void inverse(std::span<const double> y, std::span<double> x) const {
    assert(y.size() == P.size() && x.size() == P.size());
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
      : first(random_permutation(n, seed_first), false), rotation(random_blockrotation(n, seed_rotation)),
        second(random_permutation(n, seed_second), false), tmp1(n), tmp2(n) {}

  void apply(std::span<const double> x, std::span<double> y) const {
    first.apply(x, tmp1);
    rotation.apply(tmp1, tmp2);
    second.apply(tmp2, y);
  }

  void inverse(std::span<const double> y, std::span<double> x) const {
    second.inverse(y, tmp2);
    rotation.inverse(tmp2, tmp1);
    first.inverse(tmp1, x);
  }

private:
  static std::vector<size_t> random_permutation(const size_t n, const long seed) {
    std::vector<size_t> P(n);
    coco_compute_truncated_uniform_swap_permutation(P, seed, coco_get_nb_swaps(n, "bbob-largescale"),
                                                    coco_get_swap_range(n, "bbob-largescale"));
    return P;
  }

  static coco_block_diagonal_matrix_t random_blockrotation(const size_t n, const long seed) {
    size_t number_of_blocks;
    const std::vector<size_t> block_sizes = coco_get_block_sizes(&number_of_blocks, n, "bbob-largescale");
    return coco_compute_blockrotation(seed, block_sizes);
  }

  vars_permutation_t first;
  vars_blockrotation_t rotation;
  vars_permutation_t second;
  mutable std::vector<double> tmp1, tmp2; /**< @brief Intermediate values */
};

/**
 * @brief The 's_i scaling' of the Bueche-Rastrigin function: y_i = s_i x_i with s_i = sqrt(10)^(i / (n - 1)),
 * multiplied by 10 for positive x_i with even i.
 */
class vars_brs_t {
public:
  explicit vars_brs_t(const size_t n) : factors(n) {
    /* Function documentation says we should compute 10^(0.5 *
     * (i-1)/(D-1)). Instead we compute the equivalent
     * sqrt(10)^((i-1)/(D-1)) just like the legacy code.
     */
    for (size_t i = 0; i < n; ++i)
      factors[i] = pow(sqrt(10.0), static_cast<double>(i) / (static_cast<double>(n) - 1.0));
  }

  void apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == factors.size() && y.size() == factors.size());
    for (size_t i = 0; i < factors.size(); ++i) {
      double factor = factors[i];
      /* Documentation specifies odd indices and starts indexing
       * from 1, we use all even indices since C starts indexing
       * with 0.
       */
      if (x[i] > 0.0 && i % 2 == 0) {
        factor *= 10.0;
      }
      y[i] = factor * x[i];
    }
  }

  /** @brief Divides by the factor; the factor is positive, so y has the sign of x. */
  void inverse(std::span<const double> y, std::span<double> x) const {
    assert(y.size() == factors.size() && x.size() == factors.size());
    for (size_t i = 0; i < factors.size(); ++i) {
      double factor = factors[i];
      if (y[i] > 0.0 && i % 2 == 0) {
        factor *= 10.0;
      }
      x[i] = y[i] / factor;
    }
  }

private:
  std::vector<double> factors; /**< @brief sqrt(10)^(i / (n - 1)) */
};

/**
 * @brief Flips the signs of randomly chosen decision values: y_i = -x_i if the i-th of n uniform random numbers drawn
 * with the seed is below 0.5 and y_i = x_i otherwise. The transformation is its own inverse.
 */
class vars_x_hat_t {
public:
  vars_x_hat_t(const long seed, const size_t n) : flip(n) {
    std::vector<double> random_numbers(n);
    bbob2009_unif(random_numbers, seed);
    for (size_t i = 0; i < n; ++i)
      flip[i] = random_numbers[i] < 0.5;
  }

  void apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == flip.size() && y.size() == flip.size());
    for (size_t i = 0; i < flip.size(); ++i)
      y[i] = flip[i] ? -x[i] : x[i];
  }

  void inverse(std::span<const double> y, std::span<double> x) const { apply(y, x); }

private:
  std::vector<bool> flip;
};

/**
 * @brief The z^hat transformation of the Schwefel functions: y_1 = x_1 and
 * y_i = x_i + 0.25 (x_(i-1) - 2 |xopt_(i-1)|); x and y must not be the same array.
 */
class vars_z_hat_t {
public:
  explicit vars_z_hat_t(std::span<const double> xopt) : xopt(xopt.begin(), xopt.end()) {}

  void apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == xopt.size() && y.size() == xopt.size() && !xopt.empty());
    y[0] = x[0];
    for (size_t i = 1; i < xopt.size(); ++i) {
      y[i] = x[i] + 0.25 * (x[i - 1] - 2.0 * fabs(xopt[i - 1]));
    }
  }

  /** @brief Solves for x by forward substitution. */
  void inverse(std::span<const double> y, std::span<double> x) const {
    assert(y.size() == xopt.size() && x.size() == xopt.size());
    if (xopt.empty())
      return;
    x[0] = y[0];
    for (size_t i = 1; i < xopt.size(); ++i)
      x[i] = y[i] - 0.25 * (x[i - 1] - 2.0 * fabs(xopt[i - 1]));
  }

private:
  std::vector<double> xopt;
};

/**
 * @brief Rounds the decision values for the step ellipsoid: to integers where |x| > 0.5 and to multiples of
 * 1 / alpha otherwise. Not invertible.
 */
class vars_round_step_t {
public:
  explicit vars_round_step_t(const double alpha) : alpha(alpha) {}

  void apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == y.size());
    for (size_t i = 0; i < x.size(); ++i) {
      if (fabs(x[i]) > 0.5)
        y[i] = coco_double_round(x[i]);
      else
        y[i] = coco_double_round(alpha * x[i]) / alpha;
    }
  }

private:
  double alpha;
};

/**
 * @brief The penalty of the large-scale suite for x outside the region of interest [lower, upper].
 */
static double coco_boundary_penalty(std::span<const double> x, std::span<const double> lower,
                                    std::span<const double> upper) {
  assert(lower.size() == x.size() && upper.size() == x.size());
  double penalty = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
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
 * @brief The oscillating transformation T_osz of an objective value.
 */
static double coco_oscillate_objective(const double y) {
  static const double factor = 0.1;
  if (y == 0)
    return y;
  const double log_y = log(fabs(y)) / factor;
  if (y > 0)
    return pow(exp(log_y + 0.49 * (sin(log_y) + sin(0.79 * log_y))), factor);
  return -pow(exp(log_y + 0.49 * (sin(0.55 * log_y) + sin(0.31 * log_y))), factor);
}

/**
 * @brief The factor with which the large-scale suite normalizes the objective values by the dimension.
 */
static double coco_norm_by_dim_factor(const size_t n) { return bbob2009_fmin(1, 40. / static_cast<double>(n)); }
