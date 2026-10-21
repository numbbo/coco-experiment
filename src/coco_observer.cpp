/**
 * @file coco_observer.cpp
 * @brief Definitions of functions regarding COCO observers.
 */

#include "coco.h"
#include "coco_internal.h"
#include <float.h>
#include <limits.h>
#include <math.h>

namespace coco {

/**
 * @brief Number of triggers per decade of increasing evaluations (20), was: 20
 * 
 * Changing this number currently breaks the test_logger_bbob.cpp unit-test
 */
static constexpr size_t OBSERVER_DEFAULT_NUMBER_OF_EVALUATION_TRIGGERS = 20;

/**
 * @brief Number of triggers per decade of decreasing Df values (20), was: 5, then 100
 */
static constexpr size_t OBSERVER_DEFAULT_NUMBER_OF_TARGET_TRIGGERS = 20;

/**
 * @brief The type for triggers based on logarithmic target values (targets that
 * are uniformly distributed in the logarithmic space).
 *
 * The target values that trigger logging are at every
 * 10**(exponent/number_of_triggers) from positive infinity down to precision,
 * at 0, and from -precision on with step -10**(exponent/number_of_triggers)
 * until negative infinity.
 */
class ObserverLogTargets {
public:
  ObserverLogTargets() = default;
  ObserverLogTargets(const size_t number_of_triggers, const double precision);

  /** @brief Returns whether the given value triggers logging and updates the target if so. */
  bool Trigger(const double given_value);

  /** @brief Returns the currently hit target. */
  double value() const { return value_; }

private:
  int exponent_ = INT_MAX;        /**< @brief Value used to compare with the previously hit
                                    target. */
  double value_ = DBL_MAX;        /**< @brief Value of the currently hit target. */
  size_t number_of_triggers_ = 0; /**< @brief Number of target triggers between 10**i
                                    and 10**(i+1) for any i. */
  double precision_ = 0;          /**< @brief Minimal precision of interest. */
};

/**
 * @brief The type for triggers based on linear target values (targets that are
 * uniformly distributed in the linear space).
 *
 * The target values that trigger logging are at every precision * integer
 * value.
 */
class ObserverLinTargets {
public:
  ObserverLinTargets() = default;
  explicit ObserverLinTargets(const double precision) : precision_(precision) {}

  /** @brief Returns whether the given value triggers logging and updates the target if so. */
  bool Trigger(const double given_value);

  /** @brief Returns the currently hit target. */
  double value() const { return value_; }

private:
  double value_ = DBL_MAX; /**< @brief Value of the currently hit target. */
  double precision_ = 0;   /**< @brief Precision of interest. */
};

/**
 * @brief The type for triggers based on either logarithmic or linear target
 * values.
 *
 * The linear targets are always used, while the logarithmic ones are used only
 * on problems with known optima.
 */
class ObserverTargets {
public:
  ObserverTargets() = default;
  ObserverTargets(const int optima_known, const double lin_precision, const size_t number_of_targets,
                  const double log_precision);

  /** @brief Returns whether the given value triggers logging. */
  bool Trigger(const double given_value);

  /** @brief Returns the last triggered target. */
  double LastTarget() const;

private:
  bool use_log_targets_ = false;
  ObserverLinTargets lin_targets_;
  ObserverLogTargets log_targets_; /**< @brief Only used if use_log_targets is set. */
};

/**
 * @brief The type for triggers based on numbers of evaluations.
 *
 * The numbers of evaluations that trigger logging are any of the two:
 * - every 10**(exponent1/number_of_triggers) for exponent1 >= 0
 * - every base_evaluation * dimension * (10**exponent2) for exponent2 >= 0
 */
class ObserverEvaluations {
public:
  ObserverEvaluations() = default;
  ObserverEvaluations(const size_t number_evaluation_triggers, std::string_view base_evaluations,
                      const size_t dimension);

  /** @brief Returns whether the given evaluation number triggers logging. */
  bool Trigger(const size_t evaluation_number);

private:
  bool TriggerFirst(const size_t evaluation_number);
  bool TriggerSecond(const size_t evaluation_number);

  /* First trigger */
  size_t value1_ = 1;             /**< @brief The next value for the first trigger. */
  size_t exponent1_ = 0;          /**< @brief Exponent used to compute the first trigger. */
  size_t number_of_triggers_ = 0; /**< @brief Number of target triggers between 10**i
                                    and 10**(i+1) for any i. */

  /* Second trigger */
  size_t value2_ = 0;                    /**< @brief The next value for the second trigger. */
  size_t exponent2_ = 0;                 /**< @brief Exponent used to compute the second trigger. */
  std::vector<size_t> base_evaluations_; /**< @brief The base evaluation numbers used to
                                           compute the actual evaluation numbers that
                                           trigger logging. */
  size_t base_index_ = 0;                /**< @brief The next index of the base evaluations. */
  size_t dimension_ = 0;                 /**< @brief Dimension used in the calculation of the first
                                           trigger. */
};

/**
 * @brief The maximum number of evaluations to trigger logging.
 *
 * @note This is not the maximal evaluation number to be logged, but the maximal
 * number of times logging is triggered by the number of evaluations.
 */
static constexpr size_t COCO_MAX_EVALS_TO_LOG = 1000;

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
ObserverLogTargets::ObserverLogTargets(const size_t number_of_triggers, const double precision)
    : number_of_triggers_(number_of_triggers), precision_(precision) {}

/**
 * @brief Checks whether the given value should trigger logging with logarithmic
 * targets. If so, the internal values are updated.
 */
bool ObserverLogTargets::Trigger(const double given_value) {

  bool activate_trigger = false;


  double verified_value = 0;
  int current_exponent = 0;
  int adjusted_exponent = 0;

  double number_of_targets_double = static_cast<double>(number_of_triggers_);

  /* The given_value is positive or zero */
  if (given_value >= 0) {

    if (given_value == 0) {
      /* If zero, use even smaller value than precision */
      verified_value = precision_ / 10.0;
    } else if (given_value < precision_) {
      /* If close to zero, use precision instead of the given_value */
      verified_value = precision_;
    } else {
      verified_value = given_value;
    }

    current_exponent = static_cast<int>(ceil(log10(verified_value) * number_of_targets_double));

    if (current_exponent < exponent_) {
      /* Update the target information */
      exponent_ = current_exponent;
      if (given_value == 0)
        value_ = 0;
      else
        value_ = pow(10, static_cast<double>(current_exponent) / number_of_targets_double);
      activate_trigger = true;
    }
  }
  /* The given_value is negative, therefore adjustments need to be made */
  else {

    /* If close to zero, use precision instead of the given_value*/
    if (given_value > -precision_) {
      verified_value = precision_;
    } else {
      verified_value = -given_value;
    }

    /* Adjustment: use floor instead of ceil! */
    current_exponent = static_cast<int>(floor(log10(verified_value) * number_of_targets_double));

    /* Compute the adjusted_exponent in such a way, that it is always
     * diminishing in value. The adjusted exponent can only be used to verify if
     * a new target has been hit. To compute the actual target value, the
     * current_exponent needs to be used. */
    adjusted_exponent =
        2 * static_cast<int>(ceil(log10(precision_ / 10.0) * number_of_targets_double)) - current_exponent - 1;

    if (adjusted_exponent < exponent_) {
      /* Update the target information */
      exponent_ = adjusted_exponent;
      value_ = -pow(10, static_cast<double>(current_exponent) / number_of_targets_double);
      activate_trigger = true;
    }
  }

  return activate_trigger;
}

/**
 * @brief Checks whether the given value should trigger logging with linear
 * targets. If so, the internal values are updated.
 */
bool ObserverLinTargets::Trigger(const double given_value) {

  bool activate_trigger = false;

  double target_reached = DoubleRoundUpWithPrecision(given_value, precision_);
  if (target_reached < value_) {
    activate_trigger = true;
    value_ = target_reached;
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
ObserverTargets::ObserverTargets(const int optima_known, const double lin_precision, const size_t number_of_targets,
                                 const double log_precision)
    : use_log_targets_((number_of_targets > 0) && optima_known), lin_targets_(lin_precision),
      log_targets_(number_of_targets, log_precision) {}

/**
 * @brief Computes and returns whether the given value should trigger logging.
 */
bool ObserverTargets::Trigger(const double given_value) {

  bool log_trigger = use_log_targets_ && log_targets_.Trigger(given_value);
  bool lin_trigger = lin_targets_.Trigger(given_value);
  return log_trigger || lin_trigger;
}

/**
 * @brief Returns the last triggered target.
 */
double ObserverTargets::LastTarget() const {
  double last_target = lin_targets_.value();

  if (use_log_targets_) {
    double last_log_target = log_targets_.value();
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
ObserverEvaluations::ObserverEvaluations(const size_t number_evaluation_triggers, std::string_view base_evaluations,
                                         const size_t dimension)
    : number_of_triggers_(number_evaluation_triggers),
      base_evaluations_(StringParseRanges(base_evaluations, 1, 0, "base_evaluations", COCO_MAX_EVALS_TO_LOG)),
      dimension_(dimension) {

  if (this->base_evaluations_.empty())
    coco_error("ObserverEvaluations(): cannot parse base evaluations '%.*s'", static_cast<int>(base_evaluations.size()),
               base_evaluations.data());
  value2_ = dimension * this->base_evaluations_[0];
}

/**
 * @brief Computes and returns whether the given evaluation number triggers the
 * first condition of the logging based on the number of evaluations.
 *
 * The second condition is:
 * evaluation_number == 10**(exponent1/number_of_triggers)
 */
bool ObserverEvaluations::TriggerFirst(const size_t evaluation_number) {

  if (evaluation_number >= value1_) {
    /* Compute the next value for the first trigger */
    while (DoubleToSizeT(floor(pow(10, static_cast<double>(exponent1_) / static_cast<double>(number_of_triggers_)))) <=
           value1_) {
      exponent1_++;
    }
    value1_ = DoubleToSizeT(floor(pow(10, static_cast<double>(exponent1_) / static_cast<double>(number_of_triggers_))));
    return true;
  }
  return false;
}

/**
 * @brief Computes and returns whether the given evaluation number triggers the
 * second condition of the logging based on the number of evaluations.
 *
 * The second condition is:
 * evaluation_number == base_evaluation[base_index] * dimension *
 * (10**exponent2)
 */
bool ObserverEvaluations::TriggerSecond(const size_t evaluation_number) {

  if (evaluation_number >= value2_) {
    /* Compute the next value for the second trigger */
    if (base_index_ < base_evaluations_.size() - 1) {
      base_index_++;
    } else {
      base_index_ = 0;
      exponent2_++;
    }
    value2_ = DoubleToSizeT(pow(10, static_cast<double>(exponent2_)) * static_cast<double>(dimension_) *
                            static_cast<double>(base_evaluations_[base_index_]));
    return true;
  }
  return false;
}

/**
 * @brief Returns whether any of the two triggers based on the number of
 * evaluations fires.
 *
 * The numbers of evaluations that trigger logging are any of the two:
 * - every 10**(exponent1/number_of_triggers) for exponent1 >= 0
 * - every base_evaluation * dimension * (10**exponent2) for exponent2 >= 0
 */
bool ObserverEvaluations::Trigger(const size_t evaluation_number) {

  /* Both functions need to be called so that both triggers are correctly
   * updated */
  bool first = TriggerFirst(evaluation_number);
  bool second = TriggerSecond(evaluation_number);

  return first || second;
}

/**@}*/

/***********************************************************************************************************/

void Observer::SignalRestart(Problem *) {
  coco_info("coco_observer_signal_restart(): Restart signaling not supported for observer %s", observer_name.c_str());
}

} // namespace coco

#include "logger_bbob.cpp"
#include "logger_bbob_old.cpp"
#include "logger_biobj.cpp"
#include "logger_rw.cpp"
#include "logger_toy.cpp"

namespace coco {

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
static Observer *ObserverCreate(std::string_view observer_name, std::string_view observer_options) {

  Observer *observer;
  ObserverSettings observer_settings;
  std::string outer_folder = "exdata", result_folder = "default";
  int precision_x, precision_f, precision_g, log_discrete_as_int;

  size_t number_target_triggers;
  size_t number_evaluation_triggers;
  double log_target_precision, lin_target_precision;
  OptionKeys additional_option_keys; /* To be set by the chosen observer */

  /* Sets the valid keys for observer options
   * IMPORTANT: This list should be up-to-date with the code and the
   * documentation */
  OptionKeys known_option_keys = {"outer_folder",
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

  if (observer_name == "no_observer") {
    return nullptr;
  } else if (observer_name.empty()) {
    coco_warning("coco_observer(): An empty observer_name has no effect. To "
                 "prevent this warning use 'no_observer' instead");
    return nullptr;
  }

  OptionsReadString(observer_options, "outer_folder", outer_folder);
  OptionsReadString(observer_options, "result_folder", result_folder);

  /* Create the result_folder inside the outer folder */
  observer_settings.result_folder = CreateUniqueDirectory(JoinPath(outer_folder, result_folder));
  coco_info("Results will be output to folder %s", observer_settings.result_folder.c_str());

  observer_settings.observer_name = observer_name;

  observer_settings.algorithm_name = "ALG";
  OptionsReadString(observer_options, "algorithm_name", observer_settings.algorithm_name);
  OptionsReadString(observer_options, "algorithm_info", observer_settings.algorithm_info);
  OptionsReadString(observer_options, "settings", observer_settings.settings);

  number_target_triggers = OBSERVER_DEFAULT_NUMBER_OF_TARGET_TRIGGERS;
  if (OptionsReadSizeT(observer_options, "number_target_triggers", &number_target_triggers) != 0) {
    if (number_target_triggers < 0) {
      coco_warning("coco_observer(): Unsuitable observer option value "
                   "(number_target_triggers: %zu) ignored",
                   number_target_triggers);
      number_target_triggers = OBSERVER_DEFAULT_NUMBER_OF_TARGET_TRIGGERS;
    }
  }

  log_target_precision = 1e-8;
  if (OptionsReadDouble(observer_options, "log_target_precision", &log_target_precision) != 0) {
    if (log_target_precision <= 0) {
      coco_warning("coco_observer(): Unsuitable observer option value "
                   "(log_target_precision: %f) ignored",
                   log_target_precision);
      log_target_precision = 1e-8;
    }
  }

  lin_target_precision = 1e5;
  if (OptionsReadDouble(observer_options, "lin_target_precision", &lin_target_precision) != 0) {
    if (lin_target_precision <= 0) {
      coco_warning("coco_observer(): Unsuitable observer option value "
                   "(lin_target_precision: %f) ignored",
                   lin_target_precision);
      lin_target_precision = 1e5;
    }
  }

  number_evaluation_triggers = OBSERVER_DEFAULT_NUMBER_OF_EVALUATION_TRIGGERS;
  if (OptionsReadSizeT(observer_options, "number_evaluation_triggers", &number_evaluation_triggers) != 0) {
    if (number_evaluation_triggers < 4) {
      coco_warning("coco_observer(): Unsuitable observer option value "
                   "(number_evaluation_triggers: %zu) ignored",
                   number_evaluation_triggers);
      number_evaluation_triggers = OBSERVER_DEFAULT_NUMBER_OF_EVALUATION_TRIGGERS;
    }
  }

  observer_settings.base_evaluation_triggers = "1,2,5";
  OptionsReadString(observer_options, "base_evaluation_triggers", observer_settings.base_evaluation_triggers);

  precision_x = 8;
  if (OptionsReadInt(observer_options, "precision_x", &precision_x) != 0) {
    if ((precision_x < 1) || (precision_x > 32)) {
      coco_warning("coco_observer(): Unsuitable observer option value "
                   "(precision_x: %d) ignored",
                   precision_x);
      precision_x = 8;
    }
  }

  precision_f = 15;
  if (OptionsReadInt(observer_options, "precision_f", &precision_f) != 0) {
    if ((precision_f < 1) || (precision_f > 32)) {
      coco_warning("coco_observer(): Unsuitable observer option value "
                   "(precision_f: %d) ignored",
                   precision_f);
      precision_f = 15;
    }
  }

  precision_g = 3;
  if (OptionsReadInt(observer_options, "precision_g", &precision_g) != 0) {
    if ((precision_g < 1) || (precision_g > 32)) {
      coco_warning("coco_observer(): Unsuitable observer option value "
                   "(precision_g: %d) ignored",
                   precision_g);
      precision_g = 3;
    }
  }

  log_discrete_as_int = 0;
  if (OptionsReadInt(observer_options, "log_discrete_as_int", &log_discrete_as_int) != 0) {
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
             "settings: %s, number_target_triggers: %zu, log_target_precision: %e, lin_target_precision: %e, "
             "number_evaluation_triggers: %zu, base_evaluation_triggers: %s, precision_x: %d, precision_f: %d, "
             "precision_g: %d, log_discrete_as_int: %d",
             observer_settings.result_folder.c_str(), observer_settings.observer_name.c_str(),
             observer_settings.algorithm_name.c_str(), observer_settings.algorithm_info.c_str(),
             observer_settings.settings.c_str(), number_target_triggers, log_target_precision, lin_target_precision,
             number_evaluation_triggers, observer_settings.base_evaluation_triggers.c_str(), precision_x, precision_f,
             precision_g, log_discrete_as_int);

  /* Here each observer must have an entry that constructs the observer and sets the additional_option_keys */
  if (observer_name == "toy") {
    observer = new ObserverToy(observer_settings, observer_options, additional_option_keys);
  } else if (observer_name == "bbob-old") {
    observer = new ObserverBbobOld(observer_settings);
  } else if ((observer_name == "bbob") || (observer_name == "bbob-boxed") ||
             (observer_name.substr(0, 16) == "bbob-constrained") || (observer_name == "bbob-largescale") ||
             (observer_name == "bbob-mixint") || (observer_name == "bbob-noisy")) {
    observer = new ObserverBbob(observer_settings, observer_options, additional_option_keys);
  } else if ((observer_name == "bbob-biobj") || (observer_name == "bbob-biobj-ext") ||
             (observer_name == "bbob-biobj-mixint")) {
    observer = new ObserverBiobj(observer_settings, observer_options, additional_option_keys);
  } else if (observer_name == "rw") {
    observer = new ObserverRw(observer_settings, observer_options, additional_option_keys);
  } else {
    coco_warning("coco_observer(): Unknown observer %s!", observer_settings.observer_name.c_str());
    return nullptr;
  }

  /* Check for redundant option keys */
  known_option_keys.insert(known_option_keys.end(), additional_option_keys.begin(), additional_option_keys.end());
  OptionKeysWarnRedundant(observer_options, known_option_keys,
                          "coco_observer(): Some keys in observer options were ignored:\n",
                          "Valid keys for observer options are:\n");

  return observer;
}

ProblemPtr Observer::Observe(ProblemPtr problem) {

  if (problem == nullptr)
    return nullptr;

  if (!active_) {
    coco_warning("coco_problem_add_observer(): The problem will not be observed. (observer not active)");
    return problem;
  }

  return AddLogger(std::move(problem));
}

} // namespace coco

/* C interface (declared in coco.h) */

void coco_observer_free(coco_observer_t *observer) { delete FromC(observer); }

coco_observer_t *coco_observer(const char *observer_name, const char *observer_options) {
  return coco::ObserverCreate(coco::CStringView(observer_name), coco::CStringView(observer_options));
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
  if (observer == nullptr) {
    coco_warning("coco_problem_add_observer(): The problem will not be observed. (observer == NULL)");
    return FromC(problem);
  }
  return FromC(observer)->Observe(coco::ProblemPtr(FromC(problem))).release();
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

  if ((observer == nullptr) || !FromC(observer)->is_active()) {
    coco_warning("coco_problem_remove_observer(): The problem was not observed. %s",
                 observer == nullptr ? "(observer == NULL)" : "(observer not active)");
    return FromC(problem);
  }

  /* Check that we are removing the observer that is actually wrapping the
   * problem.
   *
   * This is a hack - it assumes that the name of the problem is formatted as
   * "observer_name(problem_name)". While not elegant, it does the job and is
   * better than nothing. */
  const std::string prefix = coco::RemoveFromString(FromC(problem)->name(), "(", "");
  if (prefix != FromC(observer)->name()) {
    coco_error("coco_problem_remove_observer(): trying to remove observer %s "
               "instead of %s",
               FromC(observer)->name().c_str(), prefix.c_str());
  }

  /* Keep the inner problem and free the logger */
  return ProblemTransformedFreeKeepInner(coco::ProblemPtr(FromC(problem))).release();
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
  } else if (!FromC(observer)->is_active()) {
    coco_warning("coco_observer_get_result_folder(): observer is not active, "
                 "returning empty string");
    return "";
  }
  return FromC(observer)->GetResultFolder().c_str();
}

/**
 * Invokes the logger function that stores the information about restarting the
 * algorithm (if such a function exists).
 *
 * @param problem The observed COCO problem.
 * @param observer The COCO observer that will record the restart information.
 */
void coco_observer_signal_restart(coco_observer_t *observer, coco_problem_t *problem) {

  if ((observer == nullptr) || !FromC(observer)->is_active()) {
    coco_warning("coco_observer_signal_restart(): The problem is not being observed. %s",
                 observer == nullptr ? "(observer == NULL)" : "(observer not active)");
    return;
  }

  FromC(observer)->SignalRestart(FromC(problem));
}
