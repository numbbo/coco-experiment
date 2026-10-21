/**
 * @file coco_internal.h
 * @brief Definitions of internal COCO structures and typedefs.
 *
 * These are used throughout the COCO code base but should not be used by any external code.
 */

#ifndef COCOEX_COCO_INTERNAL_H
#define COCOEX_COCO_INTERNAL_H

#include <memory>
#include <string>
#include <string_view>
#include <vector>

#ifdef __cplusplus
extern "C" {
#endif

/***********************************************************************************************************/
/**
 * @brief An owning pointer to a problem.
 *
 * Problems are created and passed around as coco_problem_ptr_t inside COCO. Ownership is only handed out as a
 * raw pointer at the C interface.
 */
using coco_problem_ptr_t = std::unique_ptr<coco_problem_t>;

/**
 * @brief The get problem function type.
 *
 * This is a template for functions that return a problem based on function, dimension and instance.
 */
using coco_get_problem_function_t = coco_problem_ptr_t (*)(const size_t function, const size_t dimension,
                                                           const size_t instance);

/**
 * @brief The option keys data type (a set of option keys used by suites and observers).
 */
using coco_option_keys_t = std::vector<std::string>;

/***********************************************************************************************************/

/**
 * @brief The COCO problem structure.
 *
 * This is one of the main structures in COCO. It contains information about a problem to be optimized.
 * Concrete problems are subclasses that override the evaluation methods. Problems can be wrapped around each
 * other (similar to the onion layers) by subclasses of coco_problem_transformed_s. Even the logger is
 * considered as just another problem wrapped around the original problem.
 *
 * The default implementations of the evaluation methods report that the method is not implemented.
 */
struct coco_problem_s {

  coco_problem_s(const size_t number_of_variables, const size_t number_of_objectives,
                 const size_t number_of_constraints);
  virtual ~coco_problem_s() = default;

  /**
   * @brief Evaluates the objective function(s) at x and stores the result in y.
   *
   * Does not check x or update any counters, see coco_evaluate_function() for that.
   */
  virtual void evaluate_function(const double *x, double *y);

  /**
   * @brief Evaluates the constraints at x and stores the result in y.
   *
   * Does not check x or update the counter of this problem, see coco_evaluate_constraint_optional_update().
   */
  virtual void evaluate_constraint(const double *x, double *y, int update_counter);

  /** @brief Evaluates the gradient of the objective function at x and stores the result in y. */
  virtual void evaluate_gradient(const double *x, double *y);

  /** @brief Logs the recommended solution x. */
  virtual void recommend_solution(const double *x);

  /**
   * @brief Sets the name and, for every variable, the same region of interest and best parameter. Assumes a
   * single-objective unconstrained problem with only continuous variables.
   */
  void init_from_scalars(std::string_view name, const double smallest, const double largest, const double best);

  /** @brief Copies the properties of the other problem, i.e. everything except its behaviour. */
  void copy_properties_from(const coco_problem_s &other);

  /** @brief Returns 1 if the best parameter is not (close to) the origin or contains NaN, 0 otherwise. */
  int best_parameter_not_zero() const;

  /**
   * @brief Returns 1 if x is finite and satisfies all constraints, 0 otherwise. Stores the constraint values
   * in constraint_values unless it is NULL.
   */
  int is_feasible(const double *x, double *constraint_values);

  /** @brief Returns the penalty of the noisy suites for x outside the region of interest. */
  double boundary_penalty(const double *x) const;

  /** @brief Sets the problem_id using formatted printing (as in printf) and checks that it is well-formed. */
  void set_id(const char *format, ...) COCO_FORMAT(2, 3);

  /** @brief Sets the problem_name using formatted printing (as in printf). */
  void set_name(const char *format, ...) COCO_FORMAT(2, 3);

  /** @brief Sets the problem_type using formatted printing (as in printf). */
  void set_type(const char *format, ...) COCO_FORMAT(2, 3);

  size_t number_of_variables;   /**< @brief Number of variables expected by the function, i.e.
                                problem dimension */
  size_t number_of_objectives;  /**< @brief Number of objectives. */
  size_t number_of_constraints; /**< @brief Number of constraints. */

  std::vector<double> smallest_values_of_interest; /**< @brief The lower bounds of the ROI in the decision space. */
  std::vector<double> largest_values_of_interest; /**< @brief The upper bounds of the ROI in the decision space. */
  size_t number_of_integer_variables;  /**< @brief Number of integer variables (if > 0, all integer variables come
                                       before any continuous ones). */

  bool is_noisy;                              /**< @brief*/
  std::vector<double> last_noise_free_values; /**< @brief Empty if the problem has no objectives. */

  std::vector<double> initial_solution; /**< @brief Initial feasible solution (empty if not defined). */
  bool is_opt_known;              /**< @brief Whether the optimal (best) value is actually known for this problem. */
  std::vector<double> best_value; /**< @brief Optimal (smallest) function value if known, otherwise a reference
                                       value defined in the problem construction. */
  std::vector<double> nadir_value; /**< @brief The nadir point (empty unless number_of_objectives > 1) */
  std::vector<double> best_parameter; /**< @brief Optimal decision vector (empty unless unique) */
  bool is_tainted; /**< @brief Is the problem tainted, i.e. has the optimal parameter setting been retrieved? */

  std::string problem_name; /**< @brief Problem name. */
  std::string problem_id;   /**< @brief Problem ID (unique in the containing suite) */
  std::string problem_type; /**< @brief Problem type (empty if not set) */

  size_t evaluations;             /**< @brief Number of objective function evaluations performed on the problem. */
  size_t evaluations_constraints; /**< @brief Number of constraint function evaluations performed on the problem. */

  /* Convenience fields for output generation */
  /* If at some point in time these arrays are changed to pointers, checks need to be added in the code to make sure
   * they are not NULL.*/

  double final_target_delta[1];       /**< @brief Final target delta. */
  double best_observed_fvalue[1];     /**< @brief The best observed value so far. */
  size_t best_observed_evaluation[1]; /**< @brief The evaluation when the best value so far was achieved. */

  /* Fields depending on the containing benchmark suite */

  coco_suite_t *suite;       /**< @brief Pointer to the containing suite (NULL if not given) */
  size_t suite_dep_index;    /**< @brief Suite-depending problem index (starting from 0) */
  size_t suite_dep_function; /**< @brief Suite-depending function */
  size_t suite_dep_instance; /**< @brief Suite-depending instance */
};

/**
 * @brief A transformed ("outer") COCO problem that wraps an inner problem.
 *
 * Transformations wrap problems like the layers of an onion. A transformed problem starts out as a copy of the
 * properties of its inner problem and by default forwards all evaluations to it. It owns the inner problem.
 */
struct coco_problem_transformed_s : public coco_problem_s {
  coco_problem_transformed_s(coco_problem_ptr_t inner, std::string_view name_prefix);

  void evaluate_function(const double *x, double *y) override;
  void evaluate_constraint(const double *x, double *y, int update_counter) override;
  void evaluate_gradient(const double *x, double *y) override;
  void recommend_solution(const double *x) override;

  coco_problem_ptr_t inner_problem; /**< @brief The inner problem, owned by this problem */
};
using coco_problem_transformed_t = coco_problem_transformed_s;

/**
 * @brief A problem constructed by stacking two problems.
 *
 * Especially useful for constructing multi-objective problems or adding constraints to a problem. It owns
 * both problems.
 */
struct coco_problem_stacked_s : public coco_problem_s {
  /**
   * @brief Stacks the two problems. Uses the given region of interest or, if it is NULL, the one of the second
   * problem.
   */
  coco_problem_stacked_s(coco_problem_ptr_t first, coco_problem_ptr_t second,
                         const double *smallest_values_of_interest = nullptr,
                         const double *largest_values_of_interest = nullptr);

  void evaluate_function(const double *x, double *y) override;
  void evaluate_constraint(const double *x, double *y, int update_counter) override;

  coco_problem_ptr_t problem1; /**< @brief The first problem (objective) */
  coco_problem_ptr_t problem2; /**< @brief The second problem (objective or constraints) */
};
using coco_problem_stacked_t = coco_problem_stacked_s;

/**
 * @brief The settings shared by all observers (read from the observer options in coco_observer()).
 */
struct coco_observer_settings_s {

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
 * problem with the corresponding logger). Every observer is a subclass that implements add_logger().
 */
struct coco_observer_s : public coco_observer_settings_s {

  explicit coco_observer_s(const coco_observer_settings_s &observer_settings)
      : coco_observer_settings_s(observer_settings), is_active(1) {}
  virtual ~coco_observer_s() = default;

  coco_observer_s(const coco_observer_s &) = delete;
  coco_observer_s &operator=(const coco_observer_s &) = delete;

  /** @brief Wraps the given problem with the logger of this observer and returns the wrapped problem. */
  virtual coco_problem_ptr_t add_logger(coco_problem_ptr_t problem) = 0;

  /**
   * @brief Signals an algorithm restart on the given (observed) problem. The default implementation reports
   * that restarts are not supported.
   */
  virtual void signal_restart(coco_problem_t *problem);

  /**
   * @brief Wraps the problem with the logger of this observer if the observer is active and returns the
   * result (the problem itself if the observer is not active, NULL for a NULL problem).
   */
  coco_problem_ptr_t observe(coco_problem_ptr_t problem);

  bool is_active; /**< @brief Whether the observer is active (the logger will log some output). */
};

/**
 * @brief The COCO suite structure.
 *
 * A suite is a collection of problems constructed by a Cartesian product of the suite's optimization
 * functions, dimensions and instances. The functions and dimensions are fixed for a suite with some name,
 * while the instances are defined dynamically. The suite can be filtered - only the chosen functions,
 * dimensions and instances will be taken into account when iterating through the suite.
 */
struct coco_suite_s {

  /**
   * @brief Sets the functions and dimensions contained in the suite, while the instances are set by
   * coco_suite_set_instance().
   */
  coco_suite_s(std::string_view suite_name, const size_t number_of_functions, const size_t number_of_dimensions,
               const size_t *dimensions, std::string_view default_instances, const int known_optima);
  /** @brief Frees the current problem (if any). */
  virtual ~coco_suite_s();

  coco_suite_s(const coco_suite_s &) = delete;
  coco_suite_s &operator=(const coco_suite_s &) = delete;

  /**
   * @brief Returns the instances of the suite for the given year. The default implementation reports an
   * error (the suite has no years defined).
   */
  virtual std::string_view get_instances_by_year(const int year) const;

  /**
   * @brief Creates the problem with the given function, dimension and instance (values, not indices).
   * get_problem_from_indices() sets the suite-dependent fields of the problem afterwards.
   */
  virtual coco_problem_ptr_t get_problem(const size_t function, const size_t dimension, const size_t instance) = 0;

  /**
   * @brief Returns the problem with the given function, dimension and instance index with its suite-dependent
   * fields set, or NULL if the problem has been filtered out.
   */
  coco_problem_ptr_t get_problem_from_indices(const size_t function_idx, const size_t dimension_idx,
                                              const size_t instance_idx);

  /** @brief Moves to the next instance that is not filtered out, returns 0 if it wrapped around. */
  int next_instance();

  /** @brief Moves to the next function that is not filtered out (resetting the instance), returns 0 if it
   * wrapped around. */
  int next_function();

  /** @brief Moves to the next dimension that is not filtered out, returns 0 if it wrapped around. */
  int next_dimension();

  std::string suite_name; /**< @brief Name of the suite. */

  std::vector<size_t> dimensions; /**< @brief The dimensions contained in the suite (0 if filtered out). */

  std::vector<size_t> functions; /**< @brief The functions contained in the suite (0 if filtered out). */

  std::string default_instances;  /**< @brief The instances contained in the suite by default. */
  std::vector<size_t> instances; /**< @brief The instances contained in the suite (0 if filtered out). */

  bool known_optima; /**< @brief Whether the problems of the suite have known optimal values */

  coco_problem_ptr_t current_problem; /**< @brief The currently tackled problem. */
  long current_dimension_idx;      /**< @brief The dimension index of the currently tackled problem. */
  long current_function_idx;       /**< @brief The function index of the currently tackled problem. */
  long current_instance_idx;       /**< @brief The instance index of the currently tackled problem. */
};

/***********************************************************************************************************/
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

/***********************************************************************************************************/

#ifdef __cplusplus
}
#endif
#endif
