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

namespace coco {

/**
 * @brief The oscillating transformation T_osz of a single value.
 */
static double ToszUv(const double xi, const double alpha) {
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
static double TasyUv(const double xi, const double beta, const size_t i, const size_t n) {
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
class VarsShift {
public:
  explicit VarsShift(std::span<const double> offset) : offset_(offset.begin(), offset.end()) {}

  void Apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == offset_.size() && y.size() == offset_.size());
    for (size_t i = 0; i < offset_.size(); ++i)
      y[i] = x[i] - offset_[i];
  }

  void Inverse(std::span<const double> y, std::span<double> x) const {
    assert(y.size() == offset_.size() && x.size() == offset_.size());
    for (size_t i = 0; i < offset_.size(); ++i)
      x[i] = y[i] + offset_[i];
  }

private:
  std::vector<double> offset_;
};

/**
 * @brief Scales the decision values: y = factor x.
 */
class VarsScale {
public:
  VarsScale(const double factor, const size_t n) : factor_(factor), n_(n) {}

  void Apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == n_ && y.size() == n_);
    for (size_t i = 0; i < n_; ++i)
      y[i] = factor_ * x[i];
  }

  /** @brief Divides by the factor. Returns false if the factor is 0. */
  bool Inverse(std::span<const double> y, std::span<double> x) const {
    assert(y.size() == n_ && x.size() == n_);
    if (factor_ == 0.)
      return false;
    for (size_t i = 0; i < n_; ++i)
      x[i] = y[i] / factor_;
    return true;
  }

private:
  double factor_;
  size_t n_;
};

/**
 * @brief Scales the decision values by alpha^(0.5 i / (n - 1)).
 */
class VarsConditioning {
public:
  VarsConditioning(const double alpha, const size_t n) : scales_(n) {
    for (size_t i = 0; i < n; ++i)
      scales_[i] = pow(alpha, 0.5 * static_cast<double>(i) / (static_cast<double>(n) - 1.0));
  }

  void Apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == scales_.size() && y.size() == scales_.size());
    for (size_t i = 0; i < scales_.size(); ++i)
      y[i] = scales_[i] * x[i];
  }

  void Inverse(std::span<const double> y, std::span<double> x) const {
    assert(y.size() == scales_.size() && x.size() == scales_.size());
    for (size_t i = 0; i < scales_.size(); ++i)
      x[i] = y[i] / scales_[i];
  }

private:
  std::vector<double> scales_; /**< @brief alpha^(0.5 i / (n - 1)) */
};

/**
 * @brief The oscillating transformation T_osz.
 */
class VarsOscillate {
public:
  void Apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == y.size());
    for (size_t i = 0; i < x.size(); ++i)
      y[i] = ToszUv(x[i], ALPHA);
  }

  /** @brief Inverts the transformation numerically. */
  void Inverse(std::span<const double> y, std::span<double> x) const {
    assert(y.size() == x.size());
    for (size_t i = 0; i < y.size(); ++i)
      x[i] = BrentInv([](const double xi) { return ToszUv(xi, ALPHA); }, y[i]);
  }

private:
  static constexpr double ALPHA = 0.1;
};

/**
 * @brief The asymmetric transformation T_asy^beta.
 */
class VarsAsymmetric {
public:
  VarsAsymmetric(const double beta, const size_t n) : beta_(beta), n_(n) {}

  void Apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == n_ && y.size() == n_);
    for (size_t i = 0; i < n_; ++i)
      y[i] = TasyUv(x[i], beta_, i, n_);
  }

  /** @brief Inverts the transformation numerically. */
  void Inverse(std::span<const double> y, std::span<double> x) const {
    assert(y.size() == n_ && x.size() == n_);
    for (size_t i = 0; i < n_; ++i)
      x[i] = BrentInv([this, i](const double xi) { return TasyUv(xi, beta_, i, n_); }, y[i]);
  }

private:
  double beta_;
  size_t n_;
};

/**
 * @brief The affine transformation y = M x + b; x and y must not be the same array.
 */
class VarsAffine {
public:
  VarsAffine(Matrix M, std::vector<double> b) : M_(std::move(M)), b_(std::move(b)) {}

  void Apply(std::span<const double> x, std::span<double> y) const { M_.MultiplyAdd(x, b_, y); }

  /** @brief Solves M x = y - b, which is M^T (y - b) if M is orthogonal. Returns false if M is singular. */
  bool Inverse(std::span<const double> y, std::span<double> x) const {
    assert(y.size() == b_.size());
    std::vector<double> shifted_y(b_.size());
    for (size_t i = 0; i < b_.size(); ++i)
      shifted_y[i] = y[i] - b_[i];
    if (M_.IsOrthogonal()) {
      M_.MultiplyTransposed(shifted_y, x);
      return true;
    }
    return M_.Solve(shifted_y, x);
  }

  const Matrix &matrix() const { return M_; }

private:
  Matrix M_; /**< @brief Has as many columns as x and as many rows as y */
  std::vector<double> b_;
};

/**
 * @brief The block rotation y = B x; x and y must not be the same array.
 */
class VarsBlockrotation {
public:
  explicit VarsBlockrotation(BlockDiagonalMatrix B) : B_(std::move(B)) {}

  void Apply(std::span<const double> x, std::span<double> y) const { B_.Multiply(x, y); }
  void Inverse(std::span<const double> y, std::span<double> x) const { B_.MultiplyTransposed(y, x); }

  const BlockDiagonalMatrix &matrix() const { return B_; }

private:
  BlockDiagonalMatrix B_;
};

/**
 * @brief Permutes the decision values: y[i] = x[P[i]], or y[P[i]] = x[i] for the inverse permutation; x and y
 * must not be the same array.
 */
class VarsPermutation {
public:
  VarsPermutation(std::span<const size_t> P, const bool is_inverse) : P_(P.begin(), P.end()), is_inverse_(is_inverse) {}

  void Apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == P_.size() && y.size() == P_.size());
    if (is_inverse_) {
      for (size_t i = 0; i < P_.size(); ++i)
        y[P_[i]] = x[i];
    } else {
      for (size_t i = 0; i < P_.size(); ++i)
        y[i] = x[P_[i]];
    }
  }

  void Inverse(std::span<const double> y, std::span<double> x) const {
    assert(y.size() == P_.size() && x.size() == P_.size());
    for (size_t i = 0; i < P_.size(); ++i) {
      if (is_inverse_)
        x[i] = y[P_[i]];
      else
        x[P_[i]] = y[i];
    }
  }

private:
  std::vector<size_t> P_; /**< @brief The permutation */
  bool is_inverse_;
};

/**
 * @brief The permuted block rotation of the large-scale suite: y = P_second B P_first x, i.e. the decision values
 * are permuted, block rotated and permuted again; x and y must not be the same array.
 *
 * The permutations are truncated uniform swaps and the block sizes those of the bbob-largescale suite.
 */
class VarsPermblockdiag {
public:
  VarsPermblockdiag(const size_t n, const long seed_first, const long seed_rotation, const long seed_second)
      : first_(RandomPermutation(n, seed_first), false), rotation_(RandomBlockrotation(n, seed_rotation)),
        second_(RandomPermutation(n, seed_second), false), tmp1_(n), tmp2_(n) {}

  void Apply(std::span<const double> x, std::span<double> y) const {
    first_.Apply(x, tmp1_);
    rotation_.Apply(tmp1_, tmp2_);
    second_.Apply(tmp2_, y);
  }

  void Inverse(std::span<const double> y, std::span<double> x) const {
    second_.Inverse(y, tmp2_);
    rotation_.Inverse(tmp2_, tmp1_);
    first_.Inverse(tmp1_, x);
  }

private:
  static std::vector<size_t> RandomPermutation(const size_t n, const long seed) {
    std::vector<size_t> P(n);
    ComputeTruncatedUniformSwapPermutation(P, seed, GetNbSwaps(n, "bbob-largescale"),
                                           GetSwapRange(n, "bbob-largescale"));
    return P;
  }

  static BlockDiagonalMatrix RandomBlockrotation(const size_t n, const long seed) {
    size_t number_of_blocks;
    const std::vector<size_t> block_sizes = GetBlockSizes(&number_of_blocks, n, "bbob-largescale");
    return ComputeBlockrotation(seed, block_sizes);
  }

  VarsPermutation first_;
  VarsBlockrotation rotation_;
  VarsPermutation second_;
  mutable std::vector<double> tmp1_, tmp2_; /**< @brief Intermediate values */
};

/**
 * @brief The 's_i scaling' of the Bueche-Rastrigin function: y_i = s_i x_i with s_i = sqrt(10)^(i / (n - 1)),
 * multiplied by 10 for positive x_i with even i.
 */
class VarsBrs {
public:
  explicit VarsBrs(const size_t n) : factors_(n) {
    /* Function documentation says we should compute 10^(0.5 *
     * (i-1)/(D-1)). Instead we compute the equivalent
     * sqrt(10)^((i-1)/(D-1)) just like the legacy code.
     */
    for (size_t i = 0; i < n; ++i)
      factors_[i] = pow(sqrt(10.0), static_cast<double>(i) / (static_cast<double>(n) - 1.0));
  }

  void Apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == factors_.size() && y.size() == factors_.size());
    for (size_t i = 0; i < factors_.size(); ++i) {
      double factor = factors_[i];
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
  void Inverse(std::span<const double> y, std::span<double> x) const {
    assert(y.size() == factors_.size() && x.size() == factors_.size());
    for (size_t i = 0; i < factors_.size(); ++i) {
      double factor = factors_[i];
      if (y[i] > 0.0 && i % 2 == 0) {
        factor *= 10.0;
      }
      x[i] = y[i] / factor;
    }
  }

private:
  std::vector<double> factors_; /**< @brief sqrt(10)^(i / (n - 1)) */
};

/**
 * @brief Flips the signs of randomly chosen decision values: y_i = -x_i if the i-th of n uniform random numbers drawn
 * with the seed is below 0.5 and y_i = x_i otherwise. The transformation is its own inverse.
 */
class VarsXHat {
public:
  VarsXHat(const long seed, const size_t n) : flip_(n) {
    std::vector<double> random_numbers(n);
    Bbob2009Unif(random_numbers, seed);
    for (size_t i = 0; i < n; ++i)
      flip_[i] = random_numbers[i] < 0.5;
  }

  void Apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == flip_.size() && y.size() == flip_.size());
    for (size_t i = 0; i < flip_.size(); ++i)
      y[i] = flip_[i] ? -x[i] : x[i];
  }

  void Inverse(std::span<const double> y, std::span<double> x) const { Apply(y, x); }

private:
  std::vector<bool> flip_;
};

/**
 * @brief The z^hat transformation of the Schwefel functions: y_1 = x_1 and
 * y_i = x_i + 0.25 (x_(i-1) - 2 |xopt_(i-1)|); x and y must not be the same array.
 */
class VarsZHat {
public:
  explicit VarsZHat(std::span<const double> xopt) : xopt_(xopt.begin(), xopt.end()) {}

  void Apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == xopt_.size() && y.size() == xopt_.size() && !xopt_.empty());
    y[0] = x[0];
    for (size_t i = 1; i < xopt_.size(); ++i) {
      y[i] = x[i] + 0.25 * (x[i - 1] - 2.0 * fabs(xopt_[i - 1]));
    }
  }

  /** @brief Solves for x by forward substitution. */
  void Inverse(std::span<const double> y, std::span<double> x) const {
    assert(y.size() == xopt_.size() && x.size() == xopt_.size());
    if (xopt_.empty())
      return;
    x[0] = y[0];
    for (size_t i = 1; i < xopt_.size(); ++i)
      x[i] = y[i] - 0.25 * (x[i - 1] - 2.0 * fabs(xopt_[i - 1]));
  }

private:
  std::vector<double> xopt_;
};

/**
 * @brief Rounds the decision values for the step ellipsoid: to integers where |x| > 0.5 and to multiples of
 * 1 / alpha otherwise. Not invertible.
 */
class VarsRoundStep {
public:
  explicit VarsRoundStep(const double alpha) : alpha_(alpha) {}

  void Apply(std::span<const double> x, std::span<double> y) const {
    assert(x.size() == y.size());
    for (size_t i = 0; i < x.size(); ++i) {
      if (fabs(x[i]) > 0.5)
        y[i] = DoubleRound(x[i]);
      else
        y[i] = DoubleRound(alpha_ * x[i]) / alpha_;
    }
  }

private:
  double alpha_;
};

/**
 * @brief The penalty of the large-scale suite for x outside the region of interest [lower, upper].
 */
static double BoundaryPenalty(std::span<const double> x, std::span<const double> lower, std::span<const double> upper) {
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
static double OscillateObjective(const double y) {
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
static double NormByDimFactor(const size_t n) { return Bbob2009Fmin(1, 40. / static_cast<double>(n)); }

} // namespace coco
