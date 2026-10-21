/**
 * @file  suite_cons_bbob_problems.cpp
 * @brief Implementation of the problems in the constrained BBOB suite.
 *
 * This suite contains 54 constrained functions in continuous domain
 * which are derived from combining 8 single-objective functions of the
 * noiseless bbob test suite with randomly-generated
 * linear constraints perturbed by nonlinear transformations.
 * Each one of the 8 functions is combined with 6 different numbers of
 * constraints: 1, 2, 6, 6+n/2, 6+n and 6+3n.
 *
 * We consider constrained optimization problems of the form
 *
 *     min f(x) s.t. g_i(x) <= 0, i=1,...,l.
 *
 * The constrained functions are built in a way such that the
 * KKT conditions are satisied at the origin, which is the initial
 * optimal solution. We then translate the constrained function
 * by a random vector in order to move the optimal solution away
 * from the origin.
 *
 * ***************************************************
 * General idea for building the constrained functions
 * ***************************************************
 *
 * 1. Choose a bbob function f whose raw version is pseudoconvex
 *    to be the objective function. (see note below)
 * 2. Remove possible nonlinear transformations from f. (see note below)
 * 3. In order to make sure that the feasible set is not empty,
 *    we choose a direction p that should be feasible.
 *    Define p as the gradient of f at the origin.
 * 4. Define the first constraint g_1(x)=a_1'*x by setting its gradient
 *    to a_1 = -p. By definition, g_1(p) < 0. Thus p is feasible
 *    for g_1.
 * 5. Generate the other constraints randomly while making sure that
 *    p remains feasible for each one.
 * 6. Apply to the whole constrained function the nonlinear transformations
 *    that were removed from the objective function in Step 2. (see note below)
 * 7. Choose a random point xopt and move the optimum away from the
 *    origin to xopt by translating the constrained function by -xopt.
 *
 * @note The pseudoconvexity in Step 1 guarantees that the KKT conditions
 *       are sufficient for optimality.
 *
 * @note The removal of possible nonlinear transformations in Step 2
 *       and posterior application in Step 6 are necessary to make
 *       sure that the KKT conditions are satisfied in the optimum
 *       - until then the origin. As explained in the documentation,
 *       the application of the nonlinear transformations in Step 6
 *       does not affect the location of the optimum.
 *
 * @note a_1 is set to -p, i.e. the negative gradient of f at the origin,
 *       in order to have the KKT conditions easily satisfied.
 *
 * @note Steps 1 and 2 are done by the ConsObjective of the problem,
 *       which leaves out the nonlinear transformations.
 *
 * @note Steps 4 and 5 are done within CLinearConsBbobConstraints().
 *
 * @note Steps 6 and 7 are done by ConsBbobProblem.
 *
 * @note The constrained Rastrigin function's construction differs a bit
 *       from the steps above. Since it is a multimodal function with
 *       well distributed local optima, we choose one of its local optima
 *       to be the global constrained optimum by adding constraints
 *       that pass through that point.
 *
 * @note The testbed provides the user an initial solution which is given
 *       by the feasible direction p scaled by some constant.
 *
 * *************************************************
 * COCO data structure for the constrained functions
 * *************************************************
 *
 * A ConsBbobProblem holds the objective function (a ConsObjective),
 * the linear constraints (see "c_linear.cpp") and the transformations
 * of steps 6 and 7, which it applies to the decision values before
 * evaluating the objective function or the constraints.
 *
 */

#include <math.h>

#include <optional>

#include "coco.h"
#include "coco_utilities.cpp"
#include "coco_random.cpp"
#include "c_linear.cpp"
#include "f_bent_cigar.cpp"
#include "f_different_powers.cpp"
#include "f_discus.cpp"
#include "f_ellipsoid.cpp"
#include "f_linear_slope.cpp"
#include "f_rastrigin.cpp"
#include "f_sphere.cpp"
#include "transform_functions.cpp"

namespace coco {

/**
 * @brief Calculates the objective function type based on the value
 *        of "function"
 */
static size_t ObjFunctionType(const size_t function) {

  size_t problems_per_obj_function_type = 6;
  return static_cast<size_t>(ceil(static_cast<double>(function) / static_cast<double>(problems_per_obj_function_type)));
}

/**
 * @brief Returns the number of linear constraints associated to the
 *        value of "function"
 */
static size_t NbOfLinearConstraints(const size_t function, const size_t dimension) {

  int problems_per_obj_function_type = 6;

  /* Map "function" value into {1, ..., problems_per_obj_function_type} */
  int p = ((static_cast<int>(function) - 1) % problems_per_obj_function_type) + 1;

  if (p == 1)
    return 1;
  else if (p == 2)
    return 2;
  else if (p == 3)
    return 6;
  else if (p == 4)
    return 6 + dimension / 2;
  else if (p == 5)
    return 6 + dimension;
  else
    return 6 + 3 * dimension;
}

/**
 * @brief Scale feasible direction depending on xopt such that
 *        xopt + feasible_direction remains in [-5, 5].
 *
 *
 */
static void FeasibleDirectionSetLength(std::span<double> feasible_direction, std::span<const double> xopt, long rseed) {
  const size_t dimension = feasible_direction.size();
  assert(xopt.size() == dimension);
  const long seed_offset = 412;    /* was sampled uniform in 0-999 */
  const double feas_shrink = 0.75; /* scale randomly between 0.75 and 1.0 */
  const double feas_bound = 5.0;

  double r[1], maxabs, maxrel;

  maxabs = maxrel = 0.0;
  for (size_t i = 0; i < dimension; ++i) {
    maxabs = DoubleMax(maxabs, fabs(xopt[i]));
    maxrel = DoubleMax(maxrel, feasible_direction[i] / (feas_bound - xopt[i]));
    maxrel = DoubleMax(maxrel, feasible_direction[i] / (-feas_bound - xopt[i]));
  }
  if (maxabs > 4.01)
    coco_warning("feasible_direction_set_length: a component of fabs(xopt) was greater than 4.01");
  if (maxabs > 5.0)
    coco_error("feasible_direction_set_length: a component of fabs(xopt) was greater than 5.0");
  Bbob2009Unif(r, rseed + seed_offset);
  VectorScale(feasible_direction, feas_shrink + r[0] * (1 - feas_shrink), /* nominator */
              maxrel);                                                    /* denominator */
}

/**
 * @brief The objective function of a constrained problem: f(x) = raw(T(x)) + fopt.
 *
 * T is the identity, a shift by xopt, a shift followed by a rotation, or (if rotate_first) a rotation followed
 * by a shift. Raw is the problem class of the raw function, e.g. Ellipsoid.
 */
template <typename Raw> class ConsObjective {
public:
  /** @brief f(x) = raw(x) + fopt */
  ConsObjective(Raw raw, const double fopt) : raw_(std::move(raw)), fopt_(fopt), z1_(raw_.number_of_variables()) {}

  /** @brief f(x) = raw(x - xopt) + fopt */
  ConsObjective(Raw raw, const double fopt, std::span<const double> xopt)
      : raw_(std::move(raw)), fopt_(fopt), shift_(xopt), z1_(raw_.number_of_variables()) {}

  /** @brief f(x) = raw(M (x - xopt)) + fopt, or f(x) = raw(M x - xopt) + fopt if rotate_first */
  ConsObjective(Raw raw, const double fopt, std::span<const double> xopt, Matrix M, const bool rotate_first)
      : raw_(std::move(raw)), fopt_(fopt), shift_(xopt),
        rotation_(VarsAffine(std::move(M), std::vector<double>(xopt.size()))), rotate_first_(rotate_first),
        z1_(raw_.number_of_variables()), z2_(raw_.number_of_variables()) {}

  double operator()(std::span<const double> x) {
    double y;
    raw_.EvaluateFunction(Transform(x), {&y, 1});
    return y + fopt_;
  }

  /** @brief Evaluates the gradient of f at x. */
  void Gradient(std::span<const double> x, std::span<double> y) {
    raw_.EvaluateGradient(Transform(x), y);
    if (!rotation_)
      return;
    /* grad_(f o g)(x), where g(x) = M * x + b, equals M^T * grad_f(g(x)) */
    const Matrix &M = rotation_->matrix();
    std::vector<double> gradient(y.size(), 0.0);
    for (size_t j = 0; j < y.size(); ++j) {
      for (size_t i = 0; i < y.size(); ++i) {
        gradient[j] += y[i] * M(i, j);
      }
    }
    std::copy(gradient.begin(), gradient.end(), y.begin());
  }

private:
  /** @brief Returns T(x). */
  std::span<const double> Transform(std::span<const double> x) {
    if (!shift_)
      return x;
    if (!rotation_) {
      shift_->Apply(x, z1_);
      return z1_;
    }
    if (rotate_first_) {
      rotation_->Apply(x, z1_);
      shift_->Apply(z1_, z2_);
    } else {
      shift_->Apply(x, z1_);
      rotation_->Apply(z1_, z2_);
    }
    return z2_;
  }

  Raw raw_;
  const double fopt_;
  const std::optional<VarsShift> shift_;
  const std::optional<VarsAffine> rotation_;
  const bool rotate_first_ = false;
  std::vector<double> z1_, z2_; /**< @brief The transformed decision values */
};

/** @brief The nonlinear transformations applied to the objective function and the constraints. */
enum class ConsNonlinear { NONE, OSCILLATE, ASYMMETRIC, ASYMMETRIC_OSCILLATE };

/**
 * @brief A constrained bbob problem.
 *
 * The objective function f and the linear constraints g_i are defined such that the optimum is at the origin.
 * The problem evaluates scale * f(T(x)) and g_i(T(x)), where T(x) = A(O(x - xopt)) with the optional
 * oscillation O and asymmetric transformation A. If not shift_objective, the objective function is evaluated as
 * scale * f(x) instead (the linear slope moves its optimum on its own).
 */
template <typename Raw> class ConsBbobProblem final : public Problem {
public:
  ConsBbobProblem(ConsObjective<Raw> objective, std::vector<LinearConstraint> constraints, std::span<const double> xopt,
                  std::span<const double> feasible_direction, const ConsNonlinear nonlinear, const bool shift_objective,
                  const double scale)
      : Problem(xopt.size(), 1, constraints.size()), objective_(std::move(objective)),
        constraints_(std::move(constraints)), shift_(xopt), asymmetric_(0.2, xopt.size()),
        oscillate_(nonlinear == ConsNonlinear::OSCILLATE || nonlinear == ConsNonlinear::ASYMMETRIC_OSCILLATE),
        asymmetric_enabled_(nonlinear == ConsNonlinear::ASYMMETRIC || nonlinear == ConsNonlinear::ASYMMETRIC_OSCILLATE),
        shift_objective_(shift_objective), scale_(scale), z_(xopt.size()) {
    for (size_t i = 0; i < number_of_variables_; ++i)
      SetRegionOfInterest(i, -5.0, 5.0);

    /* The optimum is the origin in the coordinates of the objective function and the constraints */
    const std::vector<double> origin(number_of_variables_, 0.0);
    std::vector<double> inner_x(number_of_variables_);
    NonlinearInverse(origin, inner_x);
    shift_.Inverse(inner_x, best_parameter_);
    best_value_[0] = scale_ * objective_(shift_objective_ ? std::span<const double>(origin) : best_parameter_);

    std::vector<double> initial_solution(feasible_direction.begin(), feasible_direction.end());
    if (nonlinear != ConsNonlinear::NONE)
      InverseInitialSolution(initial_solution, xopt);
    shift_.Inverse(initial_solution, initial_solution);
    SetInitialSolution(initial_solution);
  }

  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    if (shift_objective_) {
      Transform(x, z_);
      y[0] = objective_(z_);
    } else {
      y[0] = objective_(x);
    }
    y[0] *= scale_;
    AssertNotBelowBestValue(x, y);
  }

  void EvaluateConstraint(std::span<const double> x, std::span<double> y, int) override {
    Transform(x, z_);
    for (size_t j = 0; j < constraints_.size(); ++j)
      y[j] = constraints_[j](z_);
  }

private:
  /** @brief Computes z = T(x). */
  void Transform(std::span<const double> x, std::span<double> z) const {
    shift_.Apply(x, z);
    if (oscillate_)
      VarsOscillate().Apply(z, z);
    if (asymmetric_enabled_)
      asymmetric_.Apply(z, z);
  }

  /** @brief Inverts the nonlinear transformations of T. */
  void NonlinearInverse(std::span<const double> z, std::span<double> x) const {
    std::copy(z.begin(), z.end(), x.begin());
    if (asymmetric_enabled_)
      asymmetric_.Inverse(x, x);
    if (oscillate_)
      VarsOscillate().Inverse(x, x);
  }

  /**
   * @brief Maps the initial solution through the inverse of the nonlinear transformations.
   *
   * If the result plus xopt is not within (-5, 5) in every coordinate, the initial solution is halved until it
   * is. The initial solution is left unchanged if the result contains NaN.
   */
  void InverseInitialSolution(std::span<double> initial_solution, std::span<const double> xopt) const {
    const size_t n = initial_solution.size();
    const double halving_factor = .5;
    std::vector<double> halved(n), inverted(n), sol(n);
    bool is_in_bounds = true;
    double di = 0.0;

    for (size_t j = 0;; ++j) {
      for (size_t i = 0; i < n; ++i)
        halved[i] = initial_solution[i] * pow(halving_factor, static_cast<double>(j));
      NonlinearInverse(halved, inverted);
      for (size_t i = 0; i < n; ++i) {
        di = inverted[i];
        const double xi = di + xopt[i];
        is_in_bounds = (-5.0 < xi && xi < 5.0);
        if (!is_in_bounds)
          break;
        sol[i] = di;
      }
      if (is_in_bounds || IsNan(di))
        break;
    }
    if (!VectorContainsNan(sol))
      std::copy(sol.begin(), sol.end(), initial_solution.begin());
  }

  ConsObjective<Raw> objective_;
  const std::vector<LinearConstraint> constraints_;
  const VarsShift shift_;
  const VarsAsymmetric asymmetric_;
  const bool oscillate_;
  const bool asymmetric_enabled_;
  const bool shift_objective_;
  const double scale_;
  std::vector<double> z_; /**< @brief The transformed decision values */
};

/**
 * @brief Creates a constrained bbob problem from its objective function.
 *
 * The feasible direction is the gradient of the objective function at the origin, which is used to build the
 * first constraint. The objective type and the type of the constraints make up the problem type.
 */
template <typename Raw>
static ProblemPtr MakeConsBbobProblem(const size_t function, const size_t dimension, const size_t instance,
                                      const long rseed, std::span<const double> xopt, const bool active_only,
                                      const bool disguise, ConsObjective<Raw> objective, const char *objective_type,
                                      const ConsNonlinear nonlinear, const double scale,
                                      const bool shift_objective = true) {
  const std::vector<double> all_zeros(dimension, 0.0);
  std::vector<double> feasible_direction(dimension);
  objective.Gradient(all_zeros, feasible_direction);
  FeasibleDirectionSetLength(feasible_direction, xopt, rseed);

  std::vector<LinearConstraint> constraints =
      CLinearConsBbobConstraints(function, dimension, instance, NbOfLinearConstraints(function, dimension), active_only,
                                 disguise, feasible_direction);
  const char *constraints_type = constraints.size() > 1 ? "linear_linear" : "linear";

  ProblemPtr problem = std::make_unique<ConsBbobProblem<Raw>>(std::move(objective), std::move(constraints), xopt,
                                                              feasible_direction, nonlinear, shift_objective, scale);
  problem->SetType("%s_%s", objective_type, constraints_type);
  return problem;
}

/**
 * @brief Creates and returns a constrained BBOB problem.
 *
 * The scale moves the initial feasible point to a delta-f between 100 and 10000 (100 and 100,000 for the bent
 * cigar and 100 and 1000 for the different powers function).
 */
static ProblemPtr GetConsBbobProblem(std::string_view suite_name, const size_t function, const size_t dimension,
                                     const size_t instance) {
  const bool active_only = suite_name == "bbob-constrained-active-only";
  const bool disguise = suite_name != "bbob-constrained-no-disguise";

  /* Seed value used for shifting the whole constrained problem */
  std::vector<double> xopt(dimension);
  Bbob2009ComputeXopt(xopt, static_cast<long>(function + 10000 * instance));

  /* Choose a different seed value for building the objective function */
  const long rseed = static_cast<long>(function + 20000 * instance);
  const double fopt = Bbob2009ComputeFopt(function, instance);

  auto make = [&](auto objective, const char *objective_type, const ConsNonlinear nonlinear, const double scale,
                  const bool shift_objective = true) {
    return MakeConsBbobProblem(function, dimension, instance, rseed, xopt, active_only, disguise, std::move(objective),
                               objective_type, nonlinear, scale, shift_objective);
  };

  switch (ObjFunctionType(function)) {
  case 1:
    return make(ConsObjective(Sphere(dimension), fopt, BbobXopt(rseed, dimension, BbobVariant::BBOB)), "1-separable",
                ConsNonlinear::NONE, 10.);
  case 2:
    return make(ConsObjective(Ellipsoid(dimension, 1.0e6), fopt, BbobXopt(rseed, dimension, BbobVariant::BBOB)),
                "1-separable", ConsNonlinear::OSCILLATE, 1e-4);
  case 3:
    /* The linear slope is not shifted by xopt, only its constraints are */
    return make(ConsObjective(LinearSlope(dimension, BbobXopt(rseed, dimension, BbobVariant::BBOB)), fopt),
                "1-separable", ConsNonlinear::NONE, 10., false);
  case 4:
    return make(ConsObjective(Ellipsoid(dimension, 1.0e6), fopt, BbobXopt(rseed, dimension, BbobVariant::BBOB),
                              Bbob2009ComputeRotation(rseed + 1000000, dimension), false),
                "3-ill-conditioned", ConsNonlinear::OSCILLATE, 1e-4);
  case 5:
    return make(ConsObjective(Discus(dimension), fopt, BbobXopt(rseed, dimension, BbobVariant::BBOB),
                              Bbob2009ComputeRotation(rseed + 1000000, dimension), false),
                "3-ill-conditioned", ConsNonlinear::OSCILLATE, 1e-4);
  case 6:
    return make(ConsObjective(BentCigar(dimension), fopt, BbobXopt(rseed + 1000000, dimension, BbobVariant::BBOB),
                              Bbob2009ComputeRotation(rseed + 1000000, dimension), false),
                "3-ill-conditioned", ConsNonlinear::ASYMMETRIC, 1e-4);
  case 7: {
    /* The (unconstrained) optimum is set to a distance of 1e-2 to the origin. By doing so, the optimum of the
     * constrained problem is at a "reasonable" distance from the unconstrained one and, hence, the constrained
     * problem is not too easy.
     */
    std::vector<double> xopt_f = BbobXopt(rseed, dimension, BbobVariant::BBOB);
    /* CAVEAT: this implementation is not ideal for large dimensions */
    double norm = 0.0;
    for (size_t i = 0; i < dimension; ++i)
      norm += xopt_f[i] * xopt_f[i];
    norm = sqrt(norm);
    for (size_t i = 0; i < dimension; ++i)
      xopt_f[i] *= 1e-2 / norm;
    return make(ConsObjective(DifferentPowers(dimension), fopt, xopt_f,
                              Bbob2009ComputeRotation(rseed + 1000000, dimension), false),
                "3-ill-conditioned", ConsNonlinear::NONE, 1e2);
  }
  case 8: {
    std::vector<double> xopt_f(dimension);
    RastriginConsComputeXopt(xopt_f, rseed);
    return make(ConsObjective(Rastrigin(dimension), fopt, xopt_f), "1-separable", ConsNonlinear::ASYMMETRIC_OSCILLATE,
                10.);
  }
  case 9: {
    std::vector<double> xopt_f(dimension);
    RastriginConsComputeXopt(xopt_f, rseed);
    /* Rotate before shifting so that xopt is not rotated */
    return make(
        ConsObjective(Rastrigin(dimension), fopt, xopt_f, Bbob2009ComputeRotation(rseed + 1000000, dimension), true),
        "1-separable", ConsNonlinear::ASYMMETRIC_OSCILLATE, 10.);
  }
  default:
    coco_error("GetConsBbobProblem(): cannot retrieve problem f%zu instance %zu in %zuD", function, instance,
               dimension);
    return nullptr; /* Never reached */
  }
}

} // namespace coco
