#pragma once
#include "coco.h"
#include "coco_internal.h"

/**
 * @brief The type for triggers based on logarithmic target values (targets that
 * are uniformly distributed in the logarithmic space).
 *
 * The target values that trigger logging are at every
 * 10**(exponent/number_of_triggers) from positive infinity down to precision,
 * at 0, and from -precision on with step -10**(exponent/number_of_triggers)
 * until negative infinity.
 */
typedef struct coco_observer_log_targets_s {

  int exponent;              /**< @brief Value used to compare with the previously hit
                                target. */
  double value;              /**< @brief Value of the currently hit target. */
  size_t number_of_triggers; /**< @brief Number of target triggers between 10**i
                                and 10**(i+1) for any i. */
  double precision;          /**< @brief Minimal precision of interest. */

} coco_observer_log_targets_t;

/**
 * @brief The type for triggers based on linear target values (targets that are
 * uniformly distributed in the linear space).
 *
 * The target values that trigger logging are at every precision * integer
 * value.
 */
typedef struct coco_observer_lin_targets_s {

  double value;     /**< @brief Value of the currently hit target. */
  double precision; /**< @brief Precision of interest. */

} coco_observer_lin_targets_t;

/**
 * @brief The type for triggers based on either logarithmic or linear target
 * values.
 *
 * The linear targets are always used, while the logarithmic ones are used only
 * on problems with known optima.
 */
typedef struct coco_observer_targets_s {

  int use_log_targets;
  coco_observer_lin_targets_t *lin_targets;
  coco_observer_log_targets_t *log_targets;

} coco_observer_targets_t;

/**
 * @brief The type for triggers based on numbers of evaluations.
 *
 * The numbers of evaluations that trigger logging are any of the two:
 * - every 10**(exponent1/number_of_triggers) for exponent1 >= 0
 * - every base_evaluation * dimension * (10**exponent2) for exponent2 >= 0
 */
typedef struct coco_observer_evaluations_s {
  /* First trigger */
  size_t value1;             /**< @brief The next value for the first trigger. */
  size_t exponent1;          /**< @brief Exponent used to compute the first trigger. */
  size_t number_of_triggers; /**< @brief Number of target triggers between 10**i
                                and 10**(i+1) for any i. */

  /* Second trigger */
  size_t value2;            /**< @brief The next value for the second trigger. */
  size_t exponent2;         /**< @brief Exponent used to compute the second trigger. */
  size_t *base_evaluations; /**< @brief The base evaluation numbers used to
                               compute the actual evaluation numbers that
                               trigger logging. */
  size_t base_count;        /**< @brief The number of base evaluations. */
  size_t base_index;        /**< @brief The next index of the base evaluations. */
  size_t dimension;         /**< @brief Dimension used in the calculation of the first
                               trigger. */
} coco_observer_evaluations_t;

void coco_observer_signal_restart(coco_observer_t *observer, coco_problem_t *problem);
const char *coco_observer_get_result_folder(const coco_observer_t *observer);
coco_problem_t *coco_problem_remove_observer(coco_problem_t *problem, coco_observer_t *observer);
coco_problem_t *coco_problem_add_observer(coco_problem_t *problem, coco_observer_t *observer);
coco_observer_t *coco_observer(const char *observer_name, const char *observer_options);
void coco_observer_free(coco_observer_t *observer);
coco_observer_t *coco_observer_allocate(const char *result_folder, const char *observer_name,
                                        const char *algorithm_name, const char *algorithm_info,
                                        const size_t number_target_triggers, const double log_target_precision,
                                        const double lin_target_precision, const size_t number_evaluation_triggers,
                                        const char *base_evaluation_triggers, const int precision_x,
                                        const int precision_f, const int precision_g, const int log_discrete_as_int);
void coco_observer_evaluations_free(coco_observer_evaluations_t *evaluations);
int coco_observer_evaluations_trigger(coco_observer_evaluations_t *evaluations, const size_t evaluation_number);
int coco_observer_evaluations_trigger_second(coco_observer_evaluations_t *evaluations, const size_t evaluation_number);
int coco_observer_evaluations_trigger_first(coco_observer_evaluations_t *evaluations, const size_t evaluation_number);
coco_observer_evaluations_t *coco_observer_evaluations(const char *base_evaluations, const size_t dimension);
void coco_observer_targets_free(coco_observer_targets_t *targets);
double coco_observer_targets_get_last_target(coco_observer_targets_t *targets);
int coco_observer_targets_trigger(coco_observer_targets_t *targets, const double given_value);
coco_observer_targets_t *coco_observer_targets(const int optima_known, const double lin_precision,
                                               const size_t number_of_targets, const double log_precision);
int coco_observer_lin_targets_trigger(coco_observer_lin_targets_t *lin_targets, const double given_value);
coco_observer_lin_targets_t *coco_observer_lin_targets(const double precision);
int coco_observer_log_targets_trigger(coco_observer_log_targets_t *log_targets, const double given_value);
coco_observer_log_targets_t *coco_observer_log_targets(const size_t number_of_targets, const double precision);
