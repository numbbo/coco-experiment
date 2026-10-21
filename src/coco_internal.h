/**
 * @file coco_internal.h
 * @brief Definitions of internal COCO structures and typedefs.
 *
 * These are used throughout the COCO code base but should not be used by any external code.
 */

#ifndef COCOEX_COCO_INTERNAL_H
#define COCOEX_COCO_INTERNAL_H

#include <cassert>
#include <cstdio>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>

/**
 * @name The types declared in coco.h
 *
 * coco.h only declares these structures. They are empty base classes of the internal classes Problem, Suite
 * and Observer (and of Archive and RandomState, next to their definitions); the C interface converts the
 * pointers with FromC().
 */
/**@{*/
struct coco_problem_s {};
struct coco_suite_s {};
struct coco_observer_s {};
/**@}*/

namespace coco {

struct Problem;
struct Suite;
struct Observer;

/***********************************************************************************************************/
/**
 * @brief An owning pointer to a problem.
 *
 * Problems are created and passed around as ProblemPtr inside COCO. Ownership is only handed out as a
 * raw pointer at the C interface.
 */
using ProblemPtr = std::unique_ptr<Problem>;

/**
 * @brief A file that is closed when it is destroyed (or another file is opened).
 *
 * Converts to FILE * so that it can be passed to the stdio functions directly. It must not be passed to
 * fclose(), use close() instead.
 */
class File {
public:
  File() = default;
  File(const File &) = delete;
  File &operator=(const File &) = delete;
  ~File() { Close(); }

  /** @brief Opens the file at path with the given mode (see fopen()) and returns whether that worked. */
  bool Open(const std::string &path, const char *mode) {
    Close();
    file_ = fopen(path.c_str(), mode);
    return file_ != nullptr;
  }

  /** @brief Closes the file if it is open. */
  void Close() {
    if (file_ != nullptr) {
      fclose(file_);
      file_ = nullptr;
    }
  }

  operator FILE *() const { return file_; }

private:
  FILE *file_ = nullptr;
};

/**
 * @brief The get problem function type.
 *
 * This is a template for functions that return a problem based on function, dimension and instance.
 */
using GetProblemFunction = ProblemPtr (*)(const size_t function, const size_t dimension, const size_t instance);

/**
 * @brief The option keys data type (a set of option keys used by suites and observers).
 */
using OptionKeys = std::vector<std::string>;

/***********************************************************************************************************/

/**
 * @brief The COCO problem structure.
 *
 * This is one of the main structures in COCO. It contains information about a problem to be optimized.
 * Concrete problems are subclasses that override the evaluation methods. Problems can be wrapped around each
 * other (similar to the onion layers) by subclasses of TransformedProblem. Even the logger is
 * considered as just another problem wrapped around the original problem.
 *
 * The default implementations of the evaluation methods report that the method is not implemented.
 */
struct Problem : public coco_problem_s {

  Problem(const size_t number_of_variables, const size_t number_of_objectives, const size_t number_of_constraints);
  virtual ~Problem() = default;

  /**
   * @brief Evaluates the objective function(s) at x and stores the result in y.
   *
   * Does not check x or update any counters, see coco_evaluate_function() for that.
   */
  virtual void EvaluateFunction(std::span<const double> x, std::span<double> y);

  /**
   * @brief Evaluates the constraints at x and stores the result in y.
   *
   * Does not check x or update the counter of this problem, see EvaluateConstraintOptionalUpdate().
   */
  virtual void EvaluateConstraint(std::span<const double> x, std::span<double> y, int update_counter);

  /** @brief Evaluates the gradient of the objective function at x and stores the result in y. */
  virtual void EvaluateGradient(std::span<const double> x, std::span<double> y);

  /** @brief Logs the recommended solution x. */
  virtual void RecommendSolution(std::span<const double> x);

  /**
   * @brief Maps base_x from the coordinates of the base problem to the coordinates of this problem and stores
   * the result in x. Returns false if this is not possible.
   *
   * The base problem is the first problem in the chain of inner problems that is not a transformation, for a
   * problem that is not a transformation it is the problem itself (and x is a copy of base_x).
   */
  virtual bool FromBaseCoordinates(std::span<const double> base_x, std::span<double> x) const;

  /**
   * @brief Sets the name and, for every variable, the same region of interest and best parameter. Assumes a
   * single-objective unconstrained problem with only continuous variables.
   */
  void InitFromScalars(std::string_view name, const double smallest, const double largest, const double best);

  /** @brief Copies the properties of the other problem, i.e. everything except its behaviour. */
  void CopyPropertiesFrom(const Problem &other);

  /**
   * @brief Returns whether x is finite and satisfies all constraints. Stores the constraint values
   * in constraint_values unless it is empty.
   */
  bool IsFeasible(std::span<const double> x, std::span<double> constraint_values);

  /**
   * @brief Asserts that the objective value y[0] at x is not below the best value unless x is infeasible. Does
   * nothing (in particular, does not evaluate the constraints) if NDEBUG is defined.
   */
  void AssertNotBelowBestValue(std::span<const double> x, std::span<const double> y);

  /** @brief Returns the penalty of the noisy suites for x outside the region of interest. */
  double NoisyBoundaryPenalty(std::span<const double> x) const;

  /** @brief Sets the problem_id using formatted printing (as in printf) and checks that it is well-formed. */
  void SetId(const char *format, ...) COCO_FORMAT(2, 3);

  /** @brief Sets the problem_name using formatted printing (as in printf). */
  void SetName(const char *format, ...) COCO_FORMAT(2, 3);

  /** @brief Sets the problem_type using formatted printing (as in printf). */
  void SetType(const char *format, ...) COCO_FORMAT(2, 3);

  /**
   * @brief Counts an evaluation of the objective function(s) at x with the result y and keeps track of the best
   * observed feasible value.
   */
  void CountEvaluation(std::span<const double> x, std::span<const double> y);

  /** @brief Counts an evaluation of the constraints. */
  void CountConstraintEvaluation() { ++evaluations_constraints_; }

  /** @brief Resets the evaluation counters, e.g. after evaluations needed to construct the problem. */
  void ResetEvaluationCounters() {
    evaluations_ = 0;
    evaluations_constraints_ = 0;
  }

  /** @brief Marks the problem as tainted because its best parameter has been revealed. */
  void Taint() { tainted_ = true; }

  /** @brief Sets the containing suite and the suite-dependent index, function and instance. */
  void SetSuiteDependent(Suite *suite, const size_t index, const size_t function, const size_t instance);

  /** @brief Sets the suite-dependent function. */
  void set_suite_dep_function(const size_t function) { suite_dep_function_ = function; }

  /** @brief Sets whether the optimal (best) value is actually known for this problem. */
  void set_opt_known(const bool known) { opt_known_ = known; }

  /** @brief Sets the best value to the objective value(s) at the best parameter. */
  void EvaluateBestValue() { EvaluateFunction(best_parameter_, best_value_); }

  /** @brief Sets the initial solution, which should be feasible. */
  void SetInitialSolution(std::span<const double> x) {
    assert(x.size() == number_of_variables_);
    initial_solution_.assign(x.begin(), x.end());
  }

  /** @brief Sets the region of interest of the variable with the given index. */
  void SetRegionOfInterest(const size_t i, const double smallest, const double largest) {
    smallest_values_of_interest_[i] = smallest;
    largest_values_of_interest_[i] = largest;
  }

  size_t number_of_variables() const { return number_of_variables_; }
  size_t number_of_objectives() const { return number_of_objectives_; }
  size_t number_of_constraints() const { return number_of_constraints_; }
  size_t number_of_integer_variables() const { return number_of_integer_variables_; }
  const std::vector<double> &smallest_values_of_interest() const { return smallest_values_of_interest_; }
  const std::vector<double> &largest_values_of_interest() const { return largest_values_of_interest_; }
  bool is_noisy() const { return noisy_; }
  const std::vector<double> &last_noise_free_values() const { return last_noise_free_values_; }
  const std::vector<double> &initial_solution() const { return initial_solution_; }
  bool is_opt_known() const { return opt_known_; }
  const std::vector<double> &best_value() const { return best_value_; }
  const std::vector<double> &nadir_value() const { return nadir_value_; }
  /** @brief Returns the best parameter without tainting the problem (see taint()). */
  const std::vector<double> &best_parameter() const { return best_parameter_; }
  bool is_tainted() const { return tainted_; }
  const std::string &name() const { return problem_name_; }
  const std::string &id() const { return problem_id_; }
  const std::string &type() const { return problem_type_; }
  size_t evaluations() const { return evaluations_; }
  size_t evaluations_constraints() const { return evaluations_constraints_; }
  double final_target_delta() const { return final_target_delta_; }
  double best_observed_fvalue() const { return best_observed_fvalue_; }
  size_t best_observed_evaluation() const { return best_observed_evaluation_; }
  Suite *suite() const { return suite_; }
  size_t suite_dep_index() const { return suite_dep_index_; }
  size_t suite_dep_function() const { return suite_dep_function_; }
  size_t suite_dep_instance() const { return suite_dep_instance_; }

protected:
  /* The definition of the problem, which the subclasses set up in their constructors */

  size_t number_of_variables_;   /**< @brief Number of variables expected by the function, i.e.
                                problem dimension */
  size_t number_of_objectives_;  /**< @brief Number of objectives. */
  size_t number_of_constraints_; /**< @brief Number of constraints. */

  std::vector<double> smallest_values_of_interest_; /**< @brief The lower bounds of the ROI in the decision space. */
  std::vector<double> largest_values_of_interest_;  /**< @brief The upper bounds of the ROI in the decision space. */
  size_t number_of_integer_variables_; /**< @brief Number of integer variables (if > 0, all integer variables come
                                       before any continuous ones). */

  bool noisy_;                                 /**< @brief Whether the problem is noisy. */
  std::vector<double> last_noise_free_values_; /**< @brief Empty if the problem has no objectives. */

  std::vector<double> initial_solution_; /**< @brief Initial feasible solution (empty if not defined). */
  bool opt_known_;                  /**< @brief Whether the optimal (best) value is actually known for this problem. */
  std::vector<double> best_value_;  /**< @brief Optimal (smallest) function value if known, otherwise a reference
                                        value defined in the problem construction. */
  std::vector<double> nadir_value_; /**< @brief The nadir point (empty unless number_of_objectives > 1) */
  std::vector<double> best_parameter_; /**< @brief Optimal decision vector (empty unless unique) */

private:
  bool tainted_; /**< @brief Is the problem tainted, i.e. has the optimal parameter setting been retrieved? */

  std::string problem_name_; /**< @brief Problem name. */
  std::string problem_id_;   /**< @brief Problem ID (unique in the containing suite) */
  std::string problem_type_; /**< @brief Problem type (empty if not set) */

  size_t evaluations_;             /**< @brief Number of objective function evaluations performed on the problem. */
  size_t evaluations_constraints_; /**< @brief Number of constraint function evaluations performed on the problem. */

  /* Fields for output generation */
  double final_target_delta_;       /**< @brief Final target delta. */
  double best_observed_fvalue_;     /**< @brief The best observed value so far. */
  size_t best_observed_evaluation_; /**< @brief The evaluation when the best value so far was achieved. */

  /* Fields depending on the containing benchmark suite */

  Suite *suite_;              /**< @brief Pointer to the containing suite (nullptr if not given) */
  size_t suite_dep_index_;    /**< @brief Suite-depending problem index (starting from 0) */
  size_t suite_dep_function_; /**< @brief Suite-depending function */
  size_t suite_dep_instance_; /**< @brief Suite-depending instance */
};

/**
 * @brief A transformed ("outer") COCO problem that wraps an inner problem.
 *
 * Transformations wrap problems like the layers of an onion. A transformed problem starts out as a copy of the
 * properties of its inner problem and by default forwards all evaluations to it. It owns the inner problem.
 */
struct TransformedProblem : public Problem {
  TransformedProblem(ProblemPtr inner, std::string_view name_prefix);

  void EvaluateFunction(std::span<const double> x, std::span<double> y) override;
  void EvaluateConstraint(std::span<const double> x, std::span<double> y, int update_counter) override;
  void EvaluateGradient(std::span<const double> x, std::span<double> y) override;
  void RecommendSolution(std::span<const double> x) override;
  bool FromBaseCoordinates(std::span<const double> base_x, std::span<double> x) const override;

  /**
   * @brief Maps inner_x from the coordinates of the inner problem to the coordinates of this problem, i.e. inverts
   * the transformation of the decision values, and stores the result in x. Returns false if the transformation
   * is not invertible.
   */
  virtual bool Inverse(std::span<const double> inner_x, std::span<double> x) const = 0;

  Problem *inner_problem() const { return inner_problem_.get(); }

  /** @brief Takes the inner problem out of this problem, which is left without one. */
  ProblemPtr ReleaseInnerProblem() { return std::move(inner_problem_); }

protected:
  /** @brief The inverse of a transformation that does not change the decision values. */
  bool IdentityInverse(std::span<const double> inner_x, std::span<double> x) const;

  /**
   * @brief Sets the best parameter to the inverse of the best parameter of the inner problem, or to NaN if the
   * transformation is not invertible. The transformations of the decision values call it at the end of their
   * constructors.
   */
  void UpdateBestParameter();

  ProblemPtr inner_problem_; /**< @brief The inner problem, owned by this problem */
};

/**
 * @brief A problem constructed by stacking two problems.
 *
 * Especially useful for constructing multi-objective problems or adding constraints to a problem. It owns
 * both problems.
 */
struct StackedProblem : public Problem {
  /**
   * @brief Stacks the two problems. Uses the given region of interest or, if it is empty, the one of the second
   * problem.
   */
  StackedProblem(ProblemPtr first, ProblemPtr second, std::span<const double> smallest_values_of_interest = {},
                 std::span<const double> largest_values_of_interest = {});

  void EvaluateFunction(std::span<const double> x, std::span<double> y) override;
  void EvaluateConstraint(std::span<const double> x, std::span<double> y, int update_counter) override;

  Problem *problem1() const { return problem1_.get(); }
  Problem *problem2() const { return problem2_.get(); }

  /** @brief Takes the first problem out of this problem, which is left without it. */
  ProblemPtr ReleaseProblem1() { return std::move(problem1_); }

  /** @brief Takes the second problem out of this problem, which is left without it. */
  ProblemPtr ReleaseProblem2() { return std::move(problem2_); }

private:
  ProblemPtr problem1_; /**< @brief The first problem (objective) */
  ProblemPtr problem2_; /**< @brief The second problem (objective or constraints) */
};

/**
 * @brief The settings shared by all observers (read from the observer options in coco_observer()).
 */
struct ObserverSettings {

  std::string observer_name;  /**< @brief Name of the observer for identification purposes. */
  std::string result_folder;  /**< @brief Name of the result folder. */
  std::string algorithm_name; /**< @brief Name of the algorithm to be used in logger output. */
  std::string algorithm_info; /**< @brief Additional information on the algorithm to be used in logger output. */
  std::string settings;       /**< @brief Additional settings to be used in logger output. */
  size_t number_target_triggers;
  /**< @brief The number of targets between each 10**i and 10**(i+1). */
  double log_target_precision; /**< @brief The minimal precision used for logarithmic targets. */
  double lin_target_precision; /**< @brief The minimal precision used for linear targets. */
  size_t number_evaluation_triggers;
  /**< @brief The number of triggers between each 10**i and 10**(i+1) evaluation number. */
  std::string base_evaluation_triggers;
  /**< @brief The "base evaluations" used to evaluations that trigger logging. */
  int precision_x;         /**< @brief Output precision for decision variables. */
  int precision_f;         /**< @brief Output precision for function values. */
  int precision_g;         /**< @brief Output precision for constraint values. */
  bool log_discrete_as_int; /**< @brief Whether to output discrete variables in int or double format. */
};

/**
 * @brief The COCO observer structure.
 *
 * An observer observes the whole benchmark process. It is independent of suites and problems. Each time a
 * new problem of the suite is being observed, the observer initializes a new logger (wraps the observed
 * problem with the corresponding logger). Every observer is a subclass that implements AddLogger(). The
 * settings are available to the subclasses and their loggers.
 */
struct Observer : public coco_observer_s, protected ObserverSettings {

  explicit Observer(const ObserverSettings &observer_settings) : ObserverSettings(observer_settings), active_(true) {}
  virtual ~Observer() = default;

  Observer(const Observer &) = delete;
  Observer &operator=(const Observer &) = delete;

  /** @brief Wraps the given problem with the logger of this observer and returns the wrapped problem. */
  virtual ProblemPtr AddLogger(ProblemPtr problem) = 0;

  /**
   * @brief Signals an algorithm restart on the given (observed) problem. The default implementation reports
   * that restarts are not supported.
   */
  virtual void SignalRestart(Problem *problem);

  /**
   * @brief Wraps the problem with the logger of this observer if the observer is active and returns the
   * result (the problem itself if the observer is not active, nullptr for a nullptr problem).
   */
  ProblemPtr Observe(ProblemPtr problem);

  /** @brief Returns whether the observer is active (the logger will log some output). */
  bool is_active() const { return active_; }

  const std::string &name() const { return observer_name; }
  const std::string &GetResultFolder() const { return result_folder; }

protected:
  /** @brief Makes the observer inactive, e.g. because its options do not ask for any output. */
  void Deactivate() { active_ = false; }

private:
  bool active_; /**< @brief Whether the observer is active (the logger will log some output). */
};

/**
 * @brief The COCO suite structure.
 *
 * A suite is a collection of problems constructed by a Cartesian product of the suite's optimization
 * functions, dimensions and instances. The functions and dimensions are fixed for a suite with some name,
 * while the instances are defined dynamically. The suite can be filtered - only the chosen functions,
 * dimensions and instances will be taken into account when iterating through the suite.
 */
struct Suite : public coco_suite_s {

  /**
   * @brief Sets the functions and dimensions contained in the suite, while the instances are set by
   * coco_suite_set_instance().
   */
  Suite(std::string_view suite_name, const size_t number_of_functions, const size_t number_of_dimensions,
        const size_t *dimensions, std::string_view default_instances, const int known_optima);
  /** @brief Frees the current problem (if any). */
  virtual ~Suite();

  Suite(const Suite &) = delete;
  Suite &operator=(const Suite &) = delete;

  /**
   * @brief Returns the instances of the suite for the given year. The default implementation reports an
   * error (the suite has no years defined).
   */
  virtual std::string_view GetInstancesByYear(const int year) const;

  /**
   * @brief Creates the problem with the given function, dimension and instance (values, not indices).
   * GetProblemFromIndices() sets the suite-dependent fields of the problem afterwards.
   */
  virtual ProblemPtr GetProblem(const size_t function, const size_t dimension, const size_t instance) = 0;

  /**
   * @brief Returns the problem with the given function, dimension and instance index with its suite-dependent
   * fields set, or nullptr if the problem has been filtered out.
   */
  ProblemPtr GetProblemFromIndices(const size_t function_idx, const size_t dimension_idx, const size_t instance_idx);

  /** @brief Sets the instances of the suite. */
  void SetInstances(const std::vector<size_t> &instance_numbers);

  /** @brief Filters out the functions whose indices (starting from 1) are not in the given indices. */
  void FilterFunctionIndices(const std::vector<size_t> &indices);

  /** @brief Filters out the instances whose indices (starting from 1) are not in the given indices. */
  void FilterInstanceIndices(const std::vector<size_t> &indices);

  /** @brief Filters out the dimensions whose indices (starting from 1) are not in the given indices. */
  void FilterDimensionIndices(const std::vector<size_t> &indices);

  /** @brief Filters out the dimensions that are not in the given dimension numbers. */
  void FilterDimensions(const std::vector<size_t> &dimension_numbers);

  /**
   * @brief Sets the current indices such that the next call of GetNextProblem() returns the first problem
   * that is not filtered out.
   */
  void StartIteration();

  /**
   * @brief Returns the next problem that is not filtered out, wrapped with the observer unless it is nullptr,
   * or nullptr if there is no next problem. The suite owns the problem until the next call.
   */
  Problem *GetNextProblem(Observer *observer);

  /** @brief Makes the suite forget its current problem, which is then owned by the caller. */
  void ForgetCurrentProblem() { (void)current_problem_.release(); }

  const std::string &name() const { return suite_name_; }
  const std::vector<size_t> &dimensions() const { return dimensions_; }
  const std::vector<size_t> &functions() const { return functions_; }
  const std::vector<size_t> &instances() const { return instances_; }
  const std::string &default_instances() const { return default_instances_; }
  bool has_known_optima() const { return known_optima_; }

  /** @brief Sets whether the problems of the suite have known optimal values, which the tests use. */
  void set_known_optima(const bool known) { known_optima_ = known; }

private:
  /** @brief Moves to the next instance that is not filtered out, returns 0 if it wrapped around. */
  bool NextInstance();

  /** @brief Moves to the next function that is not filtered out (resetting the instance), returns 0 if it
   * wrapped around. */
  bool NextFunction();

  /** @brief Moves to the next dimension that is not filtered out, returns 0 if it wrapped around. */
  bool NextDimension();

  std::string suite_name_; /**< @brief Name of the suite. */

  std::vector<size_t> dimensions_; /**< @brief The dimensions contained in the suite (0 if filtered out). */

  std::vector<size_t> functions_; /**< @brief The functions contained in the suite (0 if filtered out). */

  std::string default_instances_; /**< @brief The instances contained in the suite by default. */
  std::vector<size_t> instances_; /**< @brief The instances contained in the suite (0 if filtered out). */

  bool known_optima_; /**< @brief Whether the problems of the suite have known optimal values */

  ProblemPtr current_problem_;         /**< @brief The currently tackled problem. */
  long current_dimension_idx_;         /**< @brief The dimension index of the currently tackled problem. */
  long current_function_idx_;          /**< @brief The function index of the currently tackled problem. */
  long current_instance_idx_;          /**< @brief The instance index of the currently tackled problem. */
};

} // namespace coco

/**
 * @name Conversions of the pointers of the C interface to the internal classes
 */
/**@{*/
inline coco::Problem *FromC(coco_problem_t *problem) { return static_cast<coco::Problem *>(problem); }
inline const coco::Problem *FromC(const coco_problem_t *problem) { return static_cast<const coco::Problem *>(problem); }
inline coco::Suite *FromC(coco_suite_t *suite) { return static_cast<coco::Suite *>(suite); }
inline const coco::Suite *FromC(const coco_suite_t *suite) { return static_cast<const coco::Suite *>(suite); }
inline coco::Observer *FromC(coco_observer_t *observer) { return static_cast<coco::Observer *>(observer); }
inline const coco::Observer *FromC(const coco_observer_t *observer) {
  return static_cast<const coco::Observer *>(observer);
}
/**@}*/

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @name Functions used by language bindings that are not part of coco.h
 *
 * They are declared here so that they keep C linkage when COCO is compiled as C++.
 */
/**@{*/
coco_problem_t *coco_get_bbob_problem(const size_t function, const size_t dimension, const size_t instance);
double coco_problem_get_best_value(const coco_problem_t *problem);
double coco_problem_get_final_target_fvalue1(const coco_problem_t *problem);
void coco_suite_forget_current_problem(coco_suite_t *suite);
/**@}*/

#ifdef __cplusplus
}
#endif
#endif
