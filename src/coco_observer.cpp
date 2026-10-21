/**
 * @file coco_observer.cpp
 * @brief Definitions of functions regarding COCO observers.
 */

#include "coco.h"
#include "coco_internal.h"
#include <float.h>
#include <limits.h>
#include <math.h>

#include <memory>

/**
 * @brief Number of triggers per decade of increasing evaluations (20), was: 20
 * 
 * Changing this number currently breaks the test_logger_bbob.cpp unit-test
 */
#define OBSERVER_DEFAULT_NUMBER_OF_EVALUATION_TRIGGERS 20

/**
 * @brief Number of triggers per decade of decreasing Df values (20), was: 5, then 100
 */
#define OBSERVER_DEFAULT_NUMBER_OF_TARGET_TRIGGERS 20

/**
 * @brief The type for triggers based on logarithmic target values (targets that
 * are uniformly distributed in the logarithmic space).
 *
 * The target values that trigger logging are at every
 * 10**(exponent/number_of_triggers) from positive infinity down to precision,
 * at 0, and from -precision on with step -10**(exponent/number_of_triggers)
 * until negative infinity.
 */
struct coco_observer_log_targets_t {

  coco_observer_log_targets_t(const size_t number_of_targets, const double precision)
      : exponent(INT_MAX), value(DBL_MAX), number_of_triggers(number_of_targets), precision(precision) {}

  int exponent;              /**< @brief Value used to compare with the previously hit
                                target. */
  double value;              /**< @brief Value of the currently hit target. */
  size_t number_of_triggers; /**< @brief Number of target triggers between 10**i
                                and 10**(i+1) for any i. */
  double precision;          /**< @brief Minimal precision of interest. */
};

/**
 * @brief The type for triggers based on linear target values (targets that are
 * uniformly distributed in the linear space).
 *
 * The target values that trigger logging are at every precision * integer
 * value.
 */
struct coco_observer_lin_targets_t {

  explicit coco_observer_lin_targets_t(const double precision) : value(DBL_MAX), precision(precision) {}

  double value;     /**< @brief Value of the currently hit target. */
  double precision; /**< @brief Precision of interest. */
};

/**
 * @brief The type for triggers based on either logarithmic or linear target
 * values.
 *
 * The linear targets are always used, while the logarithmic ones are used only
 * on problems with known optima.
 */
struct coco_observer_targets_t {

  coco_observer_targets_t(const int optima_known, const double lin_precision, const size_t number_of_targets,
                          const double log_precision)
      : use_log_targets((number_of_targets > 0) && optima_known), lin_targets(lin_precision) {
    if (use_log_targets)
      log_targets.reset(new coco_observer_log_targets_t(number_of_targets, log_precision));
  }

  int use_log_targets;
  coco_observer_lin_targets_t lin_targets;
  std::unique_ptr<coco_observer_log_targets_t> log_targets; /**< @brief NULL if log targets are not used. */
};

/**
 * @brief The type for triggers based on numbers of evaluations.
 *
 * The numbers of evaluations that trigger logging are any of the two:
 * - every 10**(exponent1/number_of_triggers) for exponent1 >= 0
 * - every base_evaluation * dimension * (10**exponent2) for exponent2 >= 0
 */
struct coco_observer_evaluations_t {

  /* First trigger */
  size_t value1;             /**< @brief The next value for the first trigger. */
  size_t exponent1;          /**< @brief Exponent used to compute the first trigger. */
  size_t number_of_triggers; /**< @brief Number of target triggers between 10**i
                                and 10**(i+1) for any i. */

  /* Second trigger */
  size_t value2;            /**< @brief The next value for the second trigger. */
  size_t exponent2;         /**< @brief Exponent used to compute the second trigger. */
  std::vector<size_t> base_evaluations; /**< @brief The base evaluation numbers used to
                                           compute the actual evaluation numbers that
                                           trigger logging. */
  size_t base_index;        /**< @brief The next index of the base evaluations. */
  size_t dimension;         /**< @brief Dimension used in the calculation of the first
                               trigger. */

};

/**
 * @brief The maximum number of evaluations to trigger logging.
 *
 * @note This is not the maximal evaluation number to be logged, but the maximal
 * number of times logging is triggered by the number of evaluations.
 */
#define COCO_MAX_EVALS_TO_LOG 1000

/***********************************************************************************************************/

/**
 * @name Methods regarding triggers based on target values
 */
/**@{*/

/**
 * @brief Creates and returns a structure containing information on logarithmic
 * targets.
 *
 * @param number_of_targets The number of targets between 10**i and 10**(i+1)
 * for each i. If 0, the logarithmic targets will not be used.
 * @param precision Minimal precision of interest.
 */
static coco_observer_log_targets_t *coco_observer_log_targets(const size_t number_of_targets, const double precision) {
  return new coco_observer_log_targets_t(number_of_targets, precision);
}

/**
 * @brief Checks whether the given value should trigger logging with logarithmic
 * targets. If so, the internal values are updated.
 */
static int coco_observer_log_targets_trigger(coco_observer_log_targets_t *log_targets, const double given_value) {

  int activate_trigger = 0;

  double number_of_targets_double;

  double verified_value = 0;
  int current_exponent = 0;
  int adjusted_exponent = 0;

  if (log_targets == nullptr) {
    return 0;
  }

  number_of_targets_double = (double)(long)log_targets->number_of_triggers;

  /* The given_value is positive or zero */
  if (given_value >= 0) {

    if (given_value == 0) {
      /* If zero, use even smaller value than precision */
      verified_value = log_targets->precision / 10.0;
    } else if (given_value < log_targets->precision) {
      /* If close to zero, use precision instead of the given_value */
      verified_value = log_targets->precision;
    } else {
      verified_value = given_value;
    }

    current_exponent = (int)(ceil(log10(verified_value) * number_of_targets_double));

    if (current_exponent < log_targets->exponent) {
      /* Update the target information */
      log_targets->exponent = current_exponent;
      if (given_value == 0)
        log_targets->value = 0;
      else
        log_targets->value = pow(10, (double)current_exponent / number_of_targets_double);
      activate_trigger = 1;
    }
  }
  /* The given_value is negative, therefore adjustments need to be made */
  else {

    /* If close to zero, use precision instead of the given_value*/
    if (given_value > -log_targets->precision) {
      verified_value = log_targets->precision;
    } else {
      verified_value = -given_value;
    }

    /* Adjustment: use floor instead of ceil! */
    current_exponent = (int)(floor(log10(verified_value) * number_of_targets_double));

    /* Compute the adjusted_exponent in such a way, that it is always
     * diminishing in value. The adjusted exponent can only be used to verify if
     * a new target has been hit. To compute the actual target value, the
     * current_exponent needs to be used. */
    adjusted_exponent =
        2 * (int)(ceil(log10(log_targets->precision / 10.0) * number_of_targets_double)) - current_exponent - 1;

    if (adjusted_exponent < log_targets->exponent) {
      /* Update the target information */
      log_targets->exponent = adjusted_exponent;
      log_targets->value = -pow(10, (double)current_exponent / number_of_targets_double);
      activate_trigger = 1;
    }
  }

  return activate_trigger;
}

/**
 * @brief Checks whether the given value should trigger logging with linear
 * targets. If so, the internal values are updated.
 */
static int coco_observer_lin_targets_trigger(coco_observer_lin_targets_t *lin_targets, const double given_value) {

  int activate_trigger = 0;
  double target_reached;

  assert(lin_targets != nullptr);

  target_reached = coco_double_round_up_with_precision(given_value, lin_targets->precision);
  if (target_reached < lin_targets->value) {
    activate_trigger = 1;
    lin_targets->value = target_reached;
  }

  return activate_trigger;
}

/**
 * @brief Creates and returns a structure containing information on triggers
 * based on linear or logarithmic target values (or both).
 *
 * If number_of_targets = 0, the logarithmic targets are not used. If not, they
 * are used only if the suite has a known optimum. The linear targets are always
 * used.
 */
static coco_observer_targets_t *coco_observer_targets(const int optima_known, const double lin_precision,
                                                      const size_t number_of_targets, const double log_precision) {
  return new coco_observer_targets_t(optima_known, lin_precision, number_of_targets, log_precision);
}

/**
 * @brief Computes and returns whether the given value should trigger logging.
 */
static int coco_observer_targets_trigger(coco_observer_targets_t *targets, const double given_value) {

  int log_trigger = coco_observer_log_targets_trigger(targets->log_targets.get(), given_value);
  int lin_trigger = coco_observer_lin_targets_trigger(&targets->lin_targets, given_value);
  return log_trigger || lin_trigger;
}

/**
 * @brief Returns the last triggered target.
 */
static double coco_observer_targets_get_last_target(const coco_observer_targets_t *targets) {
  double last_target = targets->lin_targets.value;

  if (targets->use_log_targets) {
    assert(targets->log_targets);
    double last_log_target = targets->log_targets->value;
    if (last_log_target < last_target) {
      last_target = last_log_target;
    }
  }
  return last_target;
}

/**@}*/

/***********************************************************************************************************/

/**
 * @name Methods regarding triggers based on numbers of evaluations.
 */
/**@{*/

/**
 * @brief Creates and returns a structure containing information on triggers
 * based on evaluation numbers.
 *
 * The numbers of evaluations that trigger logging are any of the two:
 * - every 10**(exponent1/number_of_triggers) for exponent1 >= 0
 * - every base_evaluation * dimension * (10**exponent2) for exponent2 >= 0
 *
 *
 * @param base_evaluations Evaluation numbers formatted as a string, which are
 * used as the base to compute the second trigger. For example, if
 * base_evaluations = "1,2,5", the logger will be triggered by evaluations
 * dim*1, dim*2, dim*5, 10*dim*1, 10*dim*2, 10*dim*5, 100*dim*1, 100*dim*2,
 * 100*dim*5, ...
 */
static coco_observer_evaluations_t *coco_observer_evaluations(const size_t number_evaluation_triggers, 
                                                              const char *base_evaluations, 
                                                              const size_t dimension) {

  coco_observer_evaluations_t *evaluations = new coco_observer_evaluations_t();

  /* First trigger */
  evaluations->value1 = 1;
  evaluations->exponent1 = 0;
  evaluations->number_of_triggers = number_evaluation_triggers;

  /* Second trigger */
  evaluations->base_evaluations =
      coco_string_parse_ranges(base_evaluations, 1, 0, "base_evaluations", COCO_MAX_EVALS_TO_LOG);
  if (evaluations->base_evaluations.empty())
    coco_error("coco_observer_evaluations(): cannot parse base evaluations '%s'", base_evaluations);
  evaluations->dimension = dimension;
  evaluations->base_index = 0;
  evaluations->value2 = dimension * evaluations->base_evaluations[0];
  evaluations->exponent2 = 0;

  return evaluations;
}

/**
 * @brief Computes and returns whether the given evaluation number triggers the
 * first condition of the logging based on the number of evaluations.
 *
 * The second condition is:
 * evaluation_number == 10**(exponent1/number_of_triggers)
 */
static int coco_observer_evaluations_trigger_first(coco_observer_evaluations_t *evaluations,
                                                   const size_t evaluation_number) {

  assert(evaluations != nullptr);

  if (evaluation_number >= evaluations->value1) {
    /* Compute the next value for the first trigger */
    while (coco_double_to_size_t(floor(pow(10, (double)evaluations->exponent1 /
                                                   (double)evaluations->number_of_triggers))) <= evaluations->value1) {
      evaluations->exponent1++;
    }
    evaluations->value1 =
        coco_double_to_size_t(floor(pow(10, (double)evaluations->exponent1 / (double)evaluations->number_of_triggers)));
    return 1;
  }
  return 0;
}

/**
 * @brief Computes and returns whether the given evaluation number triggers the
 * second condition of the logging based on the number of evaluations.
 *
 * The second condition is:
 * evaluation_number == base_evaluation[base_index] * dimension *
 * (10**exponent2)
 */
static int coco_observer_evaluations_trigger_second(coco_observer_evaluations_t *evaluations,
                                                    const size_t evaluation_number) {

  assert(evaluations != nullptr);

  if (evaluation_number >= evaluations->value2) {
    /* Compute the next value for the second trigger */
    if (evaluations->base_index < evaluations->base_evaluations.size() - 1) {
      evaluations->base_index++;
    } else {
      evaluations->base_index = 0;
      evaluations->exponent2++;
    }
    evaluations->value2 =
        coco_double_to_size_t(pow(10, (double)evaluations->exponent2) * (double)(long)evaluations->dimension *
                              (double)(long)evaluations->base_evaluations[evaluations->base_index]);
    return 1;
  }
  return 0;
}

/**
 * @brief Returns 1 if any of the two triggers based on the number of
 * evaluations equal 1 and 0 otherwise.
 *
 * The numbers of evaluations that trigger logging are any of the two:
 * - every 10**(exponent1/number_of_triggers) for exponent1 >= 0
 * - every base_evaluation * dimension * (10**exponent2) for exponent2 >= 0
 */
static int coco_observer_evaluations_trigger(coco_observer_evaluations_t *evaluations, const size_t evaluation_number) {

  /* Both functions need to be called so that both triggers are correctly
   * updated */
  int first = coco_observer_evaluations_trigger_first(evaluations, evaluation_number);
  int second = coco_observer_evaluations_trigger_second(evaluations, evaluation_number);

  return (first + second > 0) ? 1 : 0;
}

/**@}*/

/***********************************************************************************************************/

void coco_observer_free(coco_observer_t *observer) {
  delete observer;
}

void coco_observer_s::signal_restart(coco_problem_t *problem) {
  (void)problem;
  coco_info("coco_observer_signal_restart(): Restart signaling not supported for observer %s", observer_name.c_str());
}

#include "logger_bbob.cpp"
#include "logger_bbob_old.cpp"
#include "logger_biobj.cpp"
#include "logger_rw.cpp"
#include "logger_toy.cpp"

/**
 * Currently, four observers are supported:
 * - "bbob" is the observer for single-objective (both noisy and noiseless)
 * problems with known optima, which creates *.info, *.dat, *.tdat and *.rdat
 * files and logs the distance to the optimum.
 * - "bbob-biobj" is the observer for bi-objective problems, which creates
 * *.info, *.dat and *.tdat files for the given indicators, as well as an
 * archive folder with *.adat files containing nondominated solutions.
 * - "rw" is an observer for single- and bi-objective real-world problems that
 * logs all information (can be configured to long only some information) and
 * produces *.txt files (not readable by post-processing).
 * - "toy" is a simple observer that logs when a target has been hit.
 *
 * @param observer_name A string containing the name of the observer. Currently
 * supported observer names are "bbob", "bbob-biobj", "toy". Strings
 * "no_observer", "" or NULL return NULL.
 * @param observer_options A string of pairs "key: value" used to pass the
 * options to the observer. Some observer options are general, while others are
 * specific to some observers. Here we list only the general options, see
 * observer_bbob, observer_biobj and observer_toy for options of the specific
 * observers.
 * - "outer_folder: NAME" determines the outer folder for the experiment. The
 * default value is "exdata" (relative to the current folder). Absolute paths
 * are also allowed.
 * - "result_folder: NAME" determines the folder within the "exdata" folder into
 * which the results will be output. If the folder with the given name already
 * exists, first NAME_001 will be tried, then NAME_002 and so on. The default
 * value is "default".
 * - "algorithm_name: NAME", where NAME is a short name of the algorithm that
 * will be used in plots (no spaces are allowed). The default value is "ALG".
 * - "algorithm_info: STRING" stores the description of the algorithm. If it
 * contains spaces, it must be surrounded by double quotes. The default value is
 * "" (no description).
 * - "settings: STRING" stores a string with settings to be used by the post-
 * processing. The default value is "" (no settings).
 * - "number_target_triggers: VALUE" defines the number of targets between each
 * 10**i and 10**(i+1) (equally spaced in the logarithmic scale) that trigger
 * logging. The default value is 100.
 * - "log_target_precision: VALUE" defines the precision used for logarithmic
 * targets (there are no targets for abs(values) < log_target_precision). The
 * default value is 1e-8.
 * - "lin_target_precision: VALUE" defines the precision used for linear
 * targets. The default value is 1e5.
 * FIXME: This needs to be properly amended in the future!
 * - "number_evaluation_triggers: VALUE" defines the number of evaluations to be
 * logged between each 10**i and 10**(i+1). The default value is 20.
 * - "base_evaluation_triggers: VALUES" defines the base evaluations used to
 * produce an additional evaluation-based logging. The numbers of evaluations
 * that trigger logging are every base_evaluation * dimension * (10**i). For
 * example, if base_evaluation_triggers = "1,2,5", the logger will be triggered
 * by evaluations dim*1, dim*2, dim*5, 10*dim*1, 10*dim*2, 10*dim*5, 100*dim*1,
 * 100*dim*2, 100*dim*5, ... The default value is "1,2,5".
 * - "precision_x: VALUE" defines the precision used when outputting variables
 * and corresponds to the number of digits to be printed after the decimal
 * point. The default value is 8.
 * - "precision_f: VALUE" defines the precision used when outputting f values
 * and corresponds to the number of digits to be printed after the decimal
 * point. The default value is 15.
 * - "precision_g: VALUE" defines the precision used when outputting constraints
 * and corresponds to the number of digits to be printed after the decimal
 * point. The default value is 3.
 * - "log_discrete_as_int: VALUE" determines whether the values of integer
 * variables (in mixed-integer problems) are logged as integers (1) or not (0 -
 * in this case they are logged as doubles). The default value is 0.
 *
 * @return The constructed observer object or NULL if observer_name equals NULL,
 * "" or "no_observer".
 */
coco_observer_t *coco_observer(const char *observer_name, const char *observer_options) {

  coco_observer_t *observer;
  coco_observer_settings_s observer_settings;
  std::string outer_folder = "exdata", result_folder = "default";
  int precision_x, precision_f, precision_g, log_discrete_as_int;

  size_t number_target_triggers;
  size_t number_evaluation_triggers;
  double log_target_precision, lin_target_precision;
  coco_option_keys_t additional_option_keys; /* To be set by the chosen observer */

  /* Sets the valid keys for observer options
   * IMPORTANT: This list should be up-to-date with the code and the
   * documentation */
  coco_option_keys_t known_option_keys = {"outer_folder",
                              "result_folder",
                              "algorithm_name",
                              "algorithm_info",
                              "settings",
                              "number_target_triggers",
                              "log_target_precision",
                              "lin_target_precision",
                              "number_evaluation_triggers",
                              "base_evaluation_triggers",
                              "precision_x",
                              "precision_f",
                              "precision_g",
                              "log_discrete_as_int"};

  if (0 == strcmp(observer_name, "no_observer")) {
    return nullptr;
  } else if (strlen(observer_name) == 0) {
    coco_warning("coco_observer(): An empty observer_name has no effect. To "
                 "prevent this warning use 'no_observer' instead");
    return nullptr;
  }

  coco_options_read_string(observer_options, "outer_folder", outer_folder);
  coco_options_read_string(observer_options, "result_folder", result_folder);

  /* Create the result_folder inside the outer folder */
  observer_settings.result_folder = coco_create_unique_directory(coco_join_path(outer_folder, result_folder));
  coco_info("Results will be output to folder %s", observer_settings.result_folder.c_str());

  observer_settings.observer_name = observer_name;

  observer_settings.algorithm_name = "ALG";
  coco_options_read_string(observer_options, "algorithm_name", observer_settings.algorithm_name);
  coco_options_read_string(observer_options, "algorithm_info", observer_settings.algorithm_info);
  coco_options_read_string(observer_options, "settings", observer_settings.settings);

  number_target_triggers = OBSERVER_DEFAULT_NUMBER_OF_TARGET_TRIGGERS;
  if (coco_options_read_size_t(observer_options, "number_target_triggers", &number_target_triggers) != 0) {
    if (number_target_triggers < 0) {
      coco_warning("coco_observer(): Unsuitable observer option value "
                   "(number_target_triggers: %lu) ignored",
                   number_target_triggers);
      number_target_triggers = OBSERVER_DEFAULT_NUMBER_OF_TARGET_TRIGGERS;
    }
  }

  log_target_precision = 1e-8;
  if (coco_options_read_double(observer_options, "log_target_precision", &log_target_precision) != 0) {
    if (log_target_precision <= 0) {
      coco_warning("coco_observer(): Unsuitable observer option value "
                   "(log_target_precision: %f) ignored",
                   log_target_precision);
      log_target_precision = 1e-8;
    }
  }

  lin_target_precision = 1e5;
  if (coco_options_read_double(observer_options, "lin_target_precision", &lin_target_precision) != 0) {
    if (lin_target_precision <= 0) {
      coco_warning("coco_observer(): Unsuitable observer option value "
                   "(lin_target_precision: %f) ignored",
                   lin_target_precision);
      lin_target_precision = 1e5;
    }
  }

  number_evaluation_triggers = OBSERVER_DEFAULT_NUMBER_OF_EVALUATION_TRIGGERS;
  if (coco_options_read_size_t(observer_options, "number_evaluation_triggers", &number_evaluation_triggers) != 0) {
    if (number_evaluation_triggers < 4) {
      coco_warning("coco_observer(): Unsuitable observer option value "
                   "(number_evaluation_triggers: %lu) ignored",
                   number_evaluation_triggers);
      number_evaluation_triggers = OBSERVER_DEFAULT_NUMBER_OF_EVALUATION_TRIGGERS;
    }
  }

  observer_settings.base_evaluation_triggers = "1,2,5";
  coco_options_read_string(observer_options, "base_evaluation_triggers", observer_settings.base_evaluation_triggers);

  precision_x = 8;
  if (coco_options_read_int(observer_options, "precision_x", &precision_x) != 0) {
    if ((precision_x < 1) || (precision_x > 32)) {
      coco_warning("coco_observer(): Unsuitable observer option value "
                   "(precision_x: %d) ignored",
                   precision_x);
      precision_x = 8;
    }
  }

  precision_f = 15;
  if (coco_options_read_int(observer_options, "precision_f", &precision_f) != 0) {
    if ((precision_f < 1) || (precision_f > 32)) {
      coco_warning("coco_observer(): Unsuitable observer option value "
                   "(precision_f: %d) ignored",
                   precision_f);
      precision_f = 15;
    }
  }

  precision_g = 3;
  if (coco_options_read_int(observer_options, "precision_g", &precision_g) != 0) {
    if ((precision_g < 1) || (precision_g > 32)) {
      coco_warning("coco_observer(): Unsuitable observer option value "
                   "(precision_g: %d) ignored",
                   precision_g);
      precision_g = 3;
    }
  }

  log_discrete_as_int = 0;
  if (coco_options_read_int(observer_options, "log_discrete_as_int", &log_discrete_as_int) != 0) {
    if ((log_discrete_as_int < 0) || (log_discrete_as_int > 1)) {
      coco_warning("coco_observer(): Unsuitable observer option value "
                   "(log_discrete_as_int: %d) ignored",
                   log_discrete_as_int);
      log_discrete_as_int = 0;
    }
  }

  observer_settings.number_target_triggers = number_target_triggers;
  observer_settings.log_target_precision = log_target_precision;
  observer_settings.lin_target_precision = lin_target_precision;
  observer_settings.number_evaluation_triggers = number_evaluation_triggers;
  observer_settings.precision_x = precision_x;
  observer_settings.precision_f = precision_f;
  observer_settings.precision_g = precision_g;
  observer_settings.log_discrete_as_int = log_discrete_as_int;
  coco_debug("coco_observer(): result_folder: %s, observer_name: %s, algorithm_name: %s, algorithm_info: %s, "
             "settings: %s, number_target_triggers: %lu, log_target_precision: %e, lin_target_precision: %e, "
             "number_evaluation_triggers: %lu, base_evaluation_triggers: %s, precision_x: %d, precision_f: %d, "
             "precision_g: %d, log_discrete_as_int: %d",
             observer_settings.result_folder.c_str(), observer_name, observer_settings.algorithm_name.c_str(),
             observer_settings.algorithm_info.c_str(), observer_settings.settings.c_str(),
             (unsigned long)number_target_triggers, log_target_precision, lin_target_precision,
             (unsigned long)number_evaluation_triggers, observer_settings.base_evaluation_triggers.c_str(),
             precision_x, precision_f, precision_g, log_discrete_as_int);

  /* Here each observer must have an entry that constructs the observer and sets the additional_option_keys */
  if (0 == strcmp(observer_name, "toy")) {
    observer = new observer_toy_t(observer_settings, observer_options, additional_option_keys);
  } else if (0 == strcmp(observer_name, "bbob-old")) {
    observer = new observer_bbob_old_t(observer_settings);
  } else if ((0 == strcmp(observer_name, "bbob")) || (0 == strcmp(observer_name, "bbob-boxed")) ||
             (0 == strncmp(observer_name, "bbob-constrained", 16)) || (0 == strcmp(observer_name, "bbob-largescale")) ||
             (0 == strcmp(observer_name, "bbob-mixint")) || (0 == strcmp(observer_name, "bbob-noisy"))) {
    observer = new observer_bbob_t(observer_settings, observer_options, additional_option_keys);
  } else if ((0 == strcmp(observer_name, "bbob-biobj")) || (0 == strcmp(observer_name, "bbob-biobj-ext")) ||
             (0 == strcmp(observer_name, "bbob-biobj-mixint"))) {
    observer = new observer_biobj_t(observer_settings, observer_options, additional_option_keys);
  } else if (0 == strcmp(observer_name, "rw")) {
    observer = new observer_rw_t(observer_settings, observer_options, additional_option_keys);
  } else {
    coco_warning("coco_observer(): Unknown observer %s!", observer_name);
    return nullptr;
  }

  /* Check for redundant option keys */
  known_option_keys.insert(known_option_keys.end(), additional_option_keys.begin(), additional_option_keys.end());
  coco_option_keys_warn_redundant(observer_options, known_option_keys,
                                  "coco_observer(): Some keys in observer options were ignored:\n",
                                  "Valid keys for observer options are:\n");

  return observer;
}

/**
 * Wraps the observer's logger around the problem if the observer is not NULL
 * and invokes the initialization of this logger.
 *
 * @param problem The given COCO problem.
 * @param observer The COCO observer, whose logger will wrap the problem.
 *
 * @return The observed problem in the form of a new COCO problem instance or
 * the same problem if the observer is NULL.
 */
coco_problem_t *coco_problem_add_observer(coco_problem_t *problem, coco_observer_t *observer) {

  if (problem == nullptr)
    return nullptr;

  if ((observer == nullptr) || (observer->is_active == 0)) {
    coco_warning("coco_problem_add_observer(): The problem will not be observed. %s",
                 observer == nullptr ? "(observer == NULL)" : "(observer not active)");
    return problem;
  }

  return observer->add_logger(problem);
}

/**
 * Frees the observer's logger and returns the inner problem.
 *
 * @param problem The observed COCO problem.
 * @param observer The COCO observer, whose logger was wrapping the problem.
 *
 * @return The unobserved problem as a pointer to the inner problem or the same
 * problem if the problem was not observed.
 */
coco_problem_t *coco_problem_remove_observer(coco_problem_t *problem, coco_observer_t *observer) {

  coco_problem_t *problem_unobserved;
  if ((observer == nullptr) || (observer->is_active == 0)) {
    coco_warning("coco_problem_remove_observer(): The problem was not observed. %s",
                 observer == nullptr ? "(observer == NULL)" : "(observer not active)");
    return problem;
  }

  /* Check that we are removing the observer that is actually wrapping the
   * problem.
   *
   * This is a hack - it assumes that the name of the problem is formatted as
   * "observer_name(problem_name)". While not elegant, it does the job and is
   * better than nothing. */
  const std::string prefix = coco_remove_from_string(problem->problem_name, "(", "");
  if (prefix != observer->observer_name) {
    coco_error("coco_problem_remove_observer(): trying to remove observer %s "
               "instead of %s",
               observer->observer_name.c_str(), prefix.c_str());
  }

  /* Keep the inner problem and free the logger */
  problem_unobserved = coco_problem_transformed_free_keep_inner(problem);

  return problem_unobserved;
}

/**
 * Get the result folder name, which is a unique folder name constructed
 * from the result_folder option.
 *
 * @param observer The COCO observer, whose logger may be wrapping a problem.
 *
 * @return The result folder name, where the logger writes its output.
 */
const char *coco_observer_get_result_folder(const coco_observer_t *observer) {
  if (observer == nullptr) {
    coco_warning("coco_observer_get_result_folder(): no observer to get "
                 "result_folder from");
    return "";
  } else if (observer->is_active == 0) {
    coco_warning("coco_observer_get_result_folder(): observer is not active, "
                 "returning empty string");
    return "";
  }
  return observer->result_folder.c_str();
}

/**
 * Invokes the logger function that stores the information about restarting the
 * algorithm (if such a function exists).
 *
 * @param problem The observed COCO problem.
 * @param observer The COCO observer that will record the restart information.
 */
void coco_observer_signal_restart(coco_observer_t *observer, coco_problem_t *problem) {

  if ((observer == nullptr) || (observer->is_active == 0)) {
    coco_warning("coco_observer_signal_restart(): The problem is not being observed. %s",
                 observer == nullptr ? "(observer == NULL)" : "(observer not active)");
    return;
  }

  observer->signal_restart(problem);
}
