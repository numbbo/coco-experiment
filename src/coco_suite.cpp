/**
 * @file coco_suite.cpp
 * @brief Definitions of functions regarding COCO suites.
 *
 * A suite is a subclass of coco_suite_s that implements get_problem() (and
 * get_instances_by_year() if it defines instances by year). When a new suite is
 * added, the function coco_suite_intialize needs to be updated.
 *
 * @see <a href="index.html">Instructions</a> on how to write new test functions
 * and combine them into test suites.
 */

#include <time.h>

#include "coco.h"
#include "coco_internal.h"
#include "coco_utilities.cpp"

#include "suite_bbob_noisy.cpp"
#include "suite_bbob.cpp"
#include "suite_bbob_mixint.cpp"
#include "suite_biobj.cpp"
#include "suite_biobj_mixint.cpp"
#include "suite_cons_bbob.cpp"
#include "suite_largescale.cpp"
#include "suite_sbox_cost.cpp"
#include "suite_toy.cpp"

/** @brief The maximum number of different instances in a suite. */
#define COCO_MAX_INSTANCES 1000

/**
 * @brief Calls the initializer of the given suite.
 *
 * @note This function needs to be updated when a new suite is added to COCO.
 */
static coco_suite_t *coco_suite_intialize(std::string_view suite_name) {

  coco_suite_t *suite = nullptr;

  if (suite_name == "toy") {
    suite = suite_toy_initialize();
  } else if (suite_name == "bbob") {
    suite = suite_bbob_initialize();
  } else if (suite_name == "bbob-boxed") {
    suite = suite_sbox_cost_initialize();
  } else if ((suite_name == "bbob-biobj") || (suite_name == "bbob-biobj-ext")) {
    suite = suite_biobj_initialize(suite_name);
  } else if (suite_name == "bbob-largescale") {
    suite = suite_largescale_initialize();
  } else if (suite_name.substr(0, 16) == "bbob-constrained") {
    suite = suite_cons_bbob_initialize(suite_name);
  } else if (suite_name == "bbob-mixint") {
    suite = suite_bbob_mixint_initialize(suite_name);
  } else if (suite_name == "bbob-biobj-mixint") {
    suite = suite_biobj_mixint_initialize();
  } else if (suite_name == "bbob-noisy") {
    suite = suite_bbob_noisy_initialize();
  } else {
    coco_error("coco_suite_intialize(): unknown problem suite");
    return nullptr;
  }

  return suite;
}

coco_problem_ptr_t coco_suite_s::get_problem_from_indices(const size_t function_idx, const size_t dimension_idx,
                                                          const size_t instance_idx) {

  const size_t function = functions[function_idx];
  const size_t dimension = dimensions[dimension_idx];
  const size_t instance = instances[instance_idx];

  if ((function == 0) || (dimension == 0) || (instance == 0)) {
    return nullptr;
  }

  coco_problem_ptr_t problem = get_problem(function, dimension, instance);
  problem->suite_dep_function = function;
  problem->suite_dep_instance = instance;
  problem->suite_dep_index = coco_suite_encode_problem_index(this, function_idx, dimension_idx, instance_idx);
  problem->suite = this;

  return problem;
}

/**
 * @note: While a suite can contain multiple problems with equal function,
 * dimension and instance, this function always returns the first problem in the
 * suite with the given function, dimension and instance values. If the given
 * values don't correspond to a problem, the function returns NULL.
 */
coco_problem_t *coco_suite_get_problem_by_function_dimension_instance(coco_suite_t *suite, const size_t function,
                                                                      const size_t dimension, const size_t instance) {

  size_t i;
  int function_idx, dimension_idx, instance_idx;
  int found;

  found = 0;
  for (i = 0; i < suite->functions.size(); i++) {
    if (suite->functions[i] == function) {
      function_idx = (int)i;
      found = 1;
      break;
    }
  }
  if (!found)
    return nullptr;

  found = 0;
  for (i = 0; i < suite->dimensions.size(); i++) {
    if (suite->dimensions[i] == dimension) {
      dimension_idx = (int)i;
      found = 1;
      break;
    }
  }
  if (!found)
    return nullptr;

  found = 0;
  for (i = 0; i < suite->instances.size(); i++) {
    if (suite->instances[i] == instance) {
      instance_idx = (int)i;
      found = 1;
      break;
    }
  }
  if (!found)
    return nullptr;

  return suite->get_problem_from_indices((size_t)function_idx, (size_t)dimension_idx, (size_t)instance_idx).release();
}

coco_suite_s::coco_suite_s(std::string_view suite_name, const size_t number_of_functions,
                           const size_t number_of_dimensions, const size_t *dimensions,
                           std::string_view default_instances, const int known_optima)
    : suite_name(suite_name), dimensions(dimensions, dimensions + number_of_dimensions), functions(number_of_functions),
      default_instances(default_instances), known_optima(known_optima),
      /* Will be set to the first valid dimension, function and instance index by coco_suite_get_next_problem() */
      current_dimension_idx(-1), current_function_idx(-1), current_instance_idx(-1) {

  size_t i;

  assert(number_of_dimensions > 0);
  assert(number_of_functions > 0);
  for (i = 0; i < number_of_functions; i++) {
    this->functions[i] = i + 1;
    if (this->suite_name == "bbob-noisy")
      this->functions[i] = this->functions[i] + 100;
  }

  assert(!this->default_instances.empty());
  /* instances is set in coco_suite_set_instance() */
}

coco_suite_s::~coco_suite_s() {
  /* Free the current problem first, its logger might still need the suite */
  current_problem.reset();
}

std::string_view coco_suite_s::get_instances_by_year(const int year) const {
  (void)year;
  coco_error("coco_suite_get_instances_by_year(): suite '%s' has no years defined", suite_name.c_str());
  return {};
}

/**
 * @brief Sets the suite instance to the given instance_numbers.
 */
static void coco_suite_set_instance(coco_suite_t *suite, const std::vector<size_t> &instance_numbers) {

  if (instance_numbers.empty()) {
    coco_error("coco_suite_set_instance(): no instance given");
    return;
  }
  suite->instances = instance_numbers;
}

/**
 * @brief Filters the given items w.r.t. the given indices (starting from 1).
 *
 * Sets items[i - 1] to 0 for every i that cannot be found in indices (this
 * function performs the conversion from user-friendly indices starting from 1
 * to C-friendly indices starting from 0).
 */
static void coco_suite_filter_indices(std::vector<size_t> &items, const std::vector<size_t> &indices) {

  size_t i;

  for (i = 1; i <= items.size(); i++) {
    if (std::find(indices.begin(), indices.end(), i) == indices.end())
      items[i - 1] = 0;
  }
}

/**
 * @brief Filters dimensions w.r.t. the given dimension_numbers.
 *
 * Sets suite->dimensions[i] to 0 for every dimension value that cannot be found
 * in dimension_numbers.
 */
static void coco_suite_filter_dimensions(coco_suite_t *suite, const std::vector<size_t> &dimension_numbers) {

  for (size_t &dimension : suite->dimensions) {
    if (std::find(dimension_numbers.begin(), dimension_numbers.end(), dimension) == dimension_numbers.end())
      dimension = 0;
  }
}

/**
 * @param suite The given suite.
 * @param function_idx The index of the function in question (starting from 0).
 *
 * @return The function number in position function_idx in the suite. If the
 * function has been filtered out through suite_options in the coco_suite
 * function, the result is 0.
 */
size_t coco_suite_get_function_from_function_index(const coco_suite_t *suite,
                                                   const size_t function_idx) {

  if (function_idx >= suite->functions.size()) {
    coco_error("coco_suite_get_function_from_function_index(): function index "
               "exceeding the number of functions in the suite");
    return 0; /* Never reached*/
  }

  return suite->functions[function_idx];
}

/**
 * @param suite The given suite.
 * @param dimension_idx The index of the dimension in question (starting from
 * 0).
 *
 * @return The dimension number in position dimension_idx in the suite. If the
 * dimension has been filtered out through suite_options in the coco_suite
 * function, the result is 0.
 */
size_t
coco_suite_get_dimension_from_dimension_index(const coco_suite_t *suite,
                                              const size_t dimension_idx) {

  if (dimension_idx >= suite->dimensions.size()) {
    coco_error("coco_suite_get_dimension_from_dimension_index(): dimensions index exceeding the number of dimensions "
               "in the suite");
    return 0; /* Never reached*/
  }

  return suite->dimensions[dimension_idx];
}

/**
 * @param suite The given suite.
 * @param instance_idx The index of the instance in question (starting from 0).
 *
 * @return The instance number in position instance_idx in the suite. If the
 * instance has been filtered out through suite_options in the coco_suite
 * function, the result is 0.
 */
size_t coco_suite_get_instance_from_instance_index(const coco_suite_t *suite,
                                                   const size_t instance_idx) {

  if (instance_idx >= suite->instances.size()) {
    coco_error("coco_suite_get_instance_from_instance_index(): instance index "
               "exceeding the number of instances in the suite");
    return 0; /* Never reached*/
  }

  return suite->instances[instance_idx];
}

/**
 * @brief Makes the suite forget its current problem.
 *
 * Afterwards the caller owns the problem returned by the last call to coco_suite_get_next_problem() and is
 * responsible for freeing it. Used by the Rust bindings.
 */
void coco_suite_forget_current_problem(coco_suite_t *suite) {
  assert(suite != nullptr);
  (void)suite->current_problem.release(); /* Now owned by the caller */
}

void coco_suite_free(coco_suite_t *suite) {
  delete suite;
}

/**
 * Note that the problem_index depends on the number of instances a suite is
 * defined with.
 *
 * @param suite The given suite.
 * @param problem_index The index of the problem to be returned.
 *
 * @return The problem of the suite defined by problem_index (NULL if this
 * problem has been filtered out from the suite).
 */
coco_problem_t *coco_suite_get_problem(coco_suite_t *suite,
                                       const size_t problem_index) {

  size_t function_idx = 0, instance_idx = 0, dimension_idx = 0;
  coco_suite_decode_problem_index(suite, problem_index, &function_idx,
                                  &dimension_idx, &instance_idx);

  return suite->get_problem_from_indices(function_idx, dimension_idx, instance_idx).release();
}

/**
 * The number of problems in the suite is computed as a product of the number of
 * instances, number of functions and number of dimensions and therefore doesn't
 * account for any filtering done through the suite_options parameter of the
 * coco_suite function.
 *
 * @param suite The given suite.
 *
 * @return The number of problems in the suite.
 */
size_t coco_suite_get_number_of_problems(const coco_suite_t *suite) {
  return (suite->instances.size() * suite->functions.size() *
          suite->dimensions.size());
}

/**
 * @brief Returns the instances read from either a "year: YEAR" or "instances:
 * NUMBERS" string.
 *
 * If both "year" and "instances" are given, the second is ignored (and a
 * warning is raised). See the coco_suite function for more information about
 * the required format.
 */
static std::vector<size_t> coco_suite_get_instance_indices(const coco_suite_t *suite, std::string_view suite_instance) {

  int year = -1;
  std::string instances;
  std::string_view year_string;
  long year_found, instances_found;
  int parse_year = 1, parse_instances = 1;
  std::vector<size_t> result;

  year_found = coco_strfind(suite_instance, "year");
  instances_found = coco_strfind(suite_instance, "instances");

  if ((year_found < 0) && (instances_found < 0))
    return result;

  if ((year_found > 0) && (instances_found > 0)) {
    if (year_found < instances_found) {
      parse_instances = 0;
      coco_warning("coco_suite_get_instance_indices(): 'instances' suite option ignored because it follows 'year'");
    } else {
      parse_year = 0;
      coco_warning("coco_suite_get_instance_indices(): 'year' suite option "
                   "ignored because it follows 'instances'");
    }
  }

  if ((year_found >= 0) && (parse_year == 1)) {
    if (coco_options_read_int(suite_instance, "year", &(year)) != 0) {
      year_string = suite->get_instances_by_year(year);
      result = coco_string_parse_ranges(year_string, 1, 0, "instances",
                                        COCO_MAX_INSTANCES);
    } else {
      coco_warning("coco_suite_get_instance_indices(): problems parsing the "
                   "'year' suite_instance option, ignored");
    }
  }

  if ((instances_found >= 0) && (parse_instances == 1)) {
    if (coco_options_read_values(suite_instance, "instances", instances) > 0) {
      result = coco_string_parse_ranges(instances, 1, 0, "instances", COCO_MAX_INSTANCES);
    } else {
      coco_warning("coco_suite_get_instance_indices(): problems parsing the "
                   "'instance' suite_instance option, ignored");
    }
  }

  return result;
}

/**
 * @brief Iterates through the items from the current_item_id position on in
 * search for the next positive item.
 *
 * If such an item is found, current_item_id points to this item and the method
 * returns 1. If such an item cannot be found, current_item_id points to the
 * first positive item and the method returns 0.
 */
static int coco_suite_is_next_item_found(const size_t *items,
                                         const size_t number_of_items,
                                         long *current_item_id) {

  if ((*current_item_id) != (long)number_of_items - 1) {
    /* Not the last item, iterate through items */
    do {
      (*current_item_id)++;
    } while (((*current_item_id) < (long)number_of_items - 1) && (items[*current_item_id] == 0));

    assert((*current_item_id) < (long)number_of_items);
    if (items[*current_item_id] != 0) {
      /* Next item is found, return true */
      return 1;
    }
  }

  /* Next item cannot be found, move to the first good item and return false */
  *current_item_id = -1;
  do {
    (*current_item_id)++;
  } while ((*current_item_id < (long)number_of_items - 1) && (items[*current_item_id] == 0));
  if (items[*current_item_id] == 0)
    coco_error("coco_suite_is_next_item_found(): the chosen suite has no valid "
               "(positive) items");
  return 0;
}

int coco_suite_s::next_instance() {
  return coco_suite_is_next_item_found(instances.data(), instances.size(), &current_instance_idx);
}

int coco_suite_s::next_function() {
  int result = coco_suite_is_next_item_found(functions.data(), functions.size(), &current_function_idx);
  if (!result) {
    /* Reset the instances */
    current_instance_idx = -1;
    next_instance();
  }
  return result;
}

int coco_suite_s::next_dimension() {
  return coco_suite_is_next_item_found(dimensions.data(), dimensions.size(), &current_dimension_idx);
}

/**
 * Currently, seven suites are supported:
 * - "bbob" contains 24 <a
 * href="http://coco.lri.fr/downloads/download15.03/bbobdocfunctions.pdf">
 * single-objective functions</a> in 6 dimensions (2, 3, 5, 10, 20, 40)
 * - "bbob-biobj" contains 55 <a
 * href="http://numbbo.github.io/coco-doc/bbob-biobj/functions">bi-objective
 * functions</a> in 6 dimensions (2, 3, 5, 10, 20, 40)
 * - "bbob-biobj-ext" as an extension of "bbob-biobj" contains 92
 * <a href="http://numbbo.github.io/coco-doc/bbob-biobj/functions">bi-objective functions</a> in 6 dimensions
 * (2, 3, 5, 10, 20, 40)
 * - "bbob-largescale" contains 24 <a href="http://coco.lri.fr/downloads/download15.03/bbobdocfunctions.pdf">
 * single-objective functions</a> in 6 large dimensions (40, 80, 160, 320, 640, 1280)
 * - "bbob-constrained" contains 54 linearly-constrained problems, which are combinations of 8 single
 * objective functions with 6 different numbers of active linear constraints (1, 2, 10, dimension/2, dimension-1,
 * dimension+1), in 6 dimensions (2, 3, 5, 10, 20, 40).
 * - "bbob-mixint" contains 24 mixed-integer single-objective functions in 6 dimensions (2, 3, 5, 10, 20, 40)
 * - "bbob-biobj-mixint" contains 92 mixed-integer bi-objective functions in 6 dimensions (2, 3, 5, 10, 20, 40)
 * - "toy" contains 6 <a href="http://coco.lri.fr/downloads/download15.03/bbobdocfunctions.pdf">
 * single-objective functions</a> in 5 dimensions (2, 3, 5, 10, 20)
 *
 * Only the suite_name parameter needs to be non-empty. The suite_instance and
 * suite_options can be "" or NULL. In this case, default values are taken
 * (default instances of a suite are those used in the last year and the suite
 * is not filtered by default).
 *
 * @param suite_name A string containing the name of the suite. Currently
 * supported suite names are "bbob", "bbob-biobj", "bbob-biobj-ext",
 * "bbob-largescale", "bbob-constrained", and "toy".
 * @param suite_instance A string used for defining the suite instances. Two
 * ways are supported:
 * - "year: YEAR", where YEAR is the year of the BBOB workshop, includes the
 * instances (to be) used in that year's workshop;
 * - "instances: VALUES", where VALUES are instance numbers from 1 on written as
 * a comma-separated list or a range m-n.
 * @param suite_options A string of pairs "key: value" used to filter the suite
 * (especially useful for parallelizing the experiments). Supported options:
 * - "dimensions: LIST", where LIST is the list of dimensions to keep in the
 * suite (range-style syntax is not allowed here),
 * - "dimension_indices: VALUES", where VALUES is a list or a range of dimension
 * indices (starting from 1) to keep in the suite, and
 * - "function_indices: VALUES", where VALUES is a list or a range of function
 * indices (starting from 1) to keep in the suite, and
 * - "instance_indices: VALUES", where VALUES is a list or a range of instance
 * indices (starting from 1) to keep in the suite.
 *
 * @return The constructed suite object.
 */
static coco_suite_t *coco_suite_create(std::string_view suite_name, std::string_view suite_instance,
                                       std::string_view suite_options) {

  coco_suite_t *suite;
  std::vector<size_t> instances;
  std::string option_string;
  long dim_found, dim_idx_found;
  int parce_dim = 1, parce_dim_idx = 1;

  /* Sets the valid keys for suite options and suite instance */
  const coco_option_keys_t known_keys_o = {"dimensions", "dimension_indices", "function_indices", "instance_indices"};
  const coco_option_keys_t known_keys_i = {"year", "instances"};

  /* Initialize the suite */
  suite = coco_suite_intialize(suite_name);

  /* Set the instance */
  if (suite_instance.empty())
    instances = coco_suite_get_instance_indices(suite, suite->default_instances);
  else {
    instances = coco_suite_get_instance_indices(suite, suite_instance);

    if (instances.empty()) {
      /* Something wrong in the suite_instance string, use default instead */
      instances = coco_suite_get_instance_indices(suite, suite->default_instances);
    }

    /* Check for redundant option keys for suite instance */
    coco_option_keys_warn_redundant(suite_instance, known_keys_i,
                                    "coco_suite(): Some keys in suite instance were ignored:\n",
                                    "Valid keys for suite instance are:\n");
  }
  coco_suite_set_instance(suite, instances);

  /* Apply filter if any given by the suite_options */
  if (!suite_options.empty()) {
    if (coco_options_read_values(suite_options, "function_indices", option_string) > 0) {
      const std::vector<size_t> indices =
          coco_string_parse_ranges(option_string, 1, suite->functions.size(), "function_indices", COCO_MAX_INSTANCES);
      if (!indices.empty())
        coco_suite_filter_indices(suite->functions, indices);
    }

    if (coco_options_read_values(suite_options, "instance_indices", option_string) > 0) {
      const std::vector<size_t> indices =
          coco_string_parse_ranges(option_string, 1, suite->instances.size(), "instance_indices", COCO_MAX_INSTANCES);
      if (!indices.empty())
        coco_suite_filter_indices(suite->instances, indices);
    }

    dim_found = coco_strfind(suite_options, "dimensions");
    dim_idx_found = coco_strfind(suite_options, "dimension_indices");

    if ((dim_found > 0) && (dim_idx_found > 0)) {
      if (dim_found < dim_idx_found) {
        parce_dim_idx = 0;
        coco_warning("coco_suite(): 'dimension_indices' suite option ignored because it follows 'dimensions'");
      } else {
        parce_dim = 0;
        coco_warning("coco_suite(): 'dimensions' suite option ignored because "
                     "it follows 'dimension_indices'");
      }
    }

    if ((dim_idx_found >= 0) && (parce_dim_idx == 1) &&
        (coco_options_read_values(suite_options, "dimension_indices", option_string) > 0)) {
      const std::vector<size_t> indices =
          coco_string_parse_ranges(option_string, 1, suite->dimensions.size(), "dimension_indices", COCO_MAX_INSTANCES);
      if (!indices.empty())
        coco_suite_filter_indices(suite->dimensions, indices);
    }

    if ((dim_found >= 0) && (parce_dim == 1) &&
        (coco_options_read_values(suite_options, "dimensions", option_string) > 0)) {
      /* Check for disallowed characters */
      for (const char c : option_string) {
        if ((c != ',') && !isdigit((unsigned char)c)) {
          coco_warning("coco_suite(): 'dimensions' suite option ignored because of disallowed characters");
          coco_suite_free(suite);
          return nullptr;
        }
      }
      const std::vector<size_t> dimensions =
          coco_string_parse_ranges(option_string, suite->dimensions[0], suite->dimensions[suite->dimensions.size() - 1],
                                   "dimensions", COCO_MAX_INSTANCES);
      if (!dimensions.empty())
        coco_suite_filter_dimensions(suite, dimensions);
    }

    /* Check for redundant option keys for suite options */
    coco_option_keys_warn_redundant(suite_options, known_keys_o,
                                    "coco_suite(): Some keys in suite options were ignored:\n",
                                    "Valid keys for suite options are:\n");
  }

  /* Check that there are enough dimensions, functions and instances left */
  int number_of_active_dimensions = 0;
  for (size_t i = 0; i < suite->dimensions.size(); ++i) {
    number_of_active_dimensions += suite->dimensions[i] != 0;
  }
  if ((number_of_active_dimensions < 1) || (suite->dimensions.size() < 1) || (suite->functions.size() < 1) ||
      (suite->instances.size() < 1)) {
    coco_warning("coco_suite(): the suite does not contain at least one dimension, function and instance");
    coco_suite_free(suite);
    return nullptr;
  }

  /* Set the starting values of the current indices in such a way, that when the
   * instance_idx is incremented, this results in a valid problem */
  suite->next_function();
  suite->next_dimension();

  return suite;
}

coco_suite_t *coco_suite(const char *suite_name, const char *suite_instance, const char *suite_options) {
  return coco_suite_create(coco_c_string_view(suite_name), coco_c_string_view(suite_instance),
                           coco_c_string_view(suite_options));
}

/**
 * Iterates through the suite first by instances, then by functions and finally
 * by dimensions. The instances/functions/dimensions that have been filtered out
 * using the suite_options of the coco_suite function are skipped. Outputs some
 * information regarding the current place in the iteration. The returned
 * problem is wrapped with the observer. If the observer is NULL, the returned
 * problem is unobserved.
 *
 * @param suite The given suite.
 * @param observer The observer used to wrap the problem. If NULL, the problem
 * is returned unobserved.
 *
 * @returns The next problem of the suite or NULL if there is no next problem
 * left.
 */
coco_problem_t *coco_suite_get_next_problem(coco_suite_t *suite, coco_observer_t *observer) {

  size_t function_idx;
  size_t dimension_idx;
  size_t instance_idx;
  coco_problem_ptr_t problem;

  long previous_function_idx;
  long previous_dimension_idx;
  long previous_instance_idx;

  assert(suite != nullptr);

  previous_function_idx = suite->current_function_idx;
  previous_dimension_idx = suite->current_dimension_idx;
  previous_instance_idx = suite->current_instance_idx;

  /* Iterate through the suite by instances, then functions and lastly dimensions in search for the next
   * problem. Note that these functions set the values of suite fields current_instance_idx,
   * current_function_idx and current_dimension_idx. */
  if (!suite->next_instance() && !suite->next_function() && !suite->next_dimension()) {
    coco_info_partial("done\n");
    return nullptr;
  }

  suite->current_problem.reset();

  assert(suite->current_function_idx >= 0);
  assert(suite->current_dimension_idx >= 0);
  assert(suite->current_instance_idx >= 0);

  function_idx = (size_t)suite->current_function_idx;
  dimension_idx = (size_t)suite->current_dimension_idx;
  instance_idx = (size_t)suite->current_instance_idx;

  problem = suite->get_problem_from_indices(function_idx, dimension_idx, instance_idx);
  if (observer != nullptr)
    problem = observer->observe(std::move(problem));
  suite->current_problem = std::move(problem);

  /* Output information regarding the current place in the iteration */
  if (coco_log_level >= COCO_INFO) {
    if (((long)dimension_idx != previous_dimension_idx) || (previous_instance_idx < 0)) {
      /* A new dimension started */
      const std::string time_string = coco_current_time_get_string();
      if (dimension_idx > 0)
        coco_info_partial("done\n");
      else
        coco_info_partial("\n");
      coco_info_partial("COCO INFO: %s, d=%lu, running: f%02lu", time_string.c_str(),
                        (unsigned long)suite->dimensions[dimension_idx], (unsigned long)suite->functions[function_idx]);
    } else if ((long)function_idx != previous_function_idx) {
      /* A new function started */
      coco_info_partial("f%02lu", (unsigned long)suite->functions[function_idx]);
    }
    /* One dot for each instance */
    coco_info_partial(".");
  }

  return suite->current_problem.get();
}

/* See coco.h for more information on encoding and decoding problem index */

/**
 * @param suite The suite.
 * @param function_idx Index of the function (starting with 0).
 * @param dimension_idx Index of the dimension (starting with 0).
 * @param instance_idx Index of the instance (starting with 0).
 *
 * @return The problem index in the suite computed from function_idx,
 * dimension_idx and instance_idx.
 */
size_t coco_suite_encode_problem_index(const coco_suite_t *suite, const size_t function_idx, const size_t dimension_idx,
                                       const size_t instance_idx) {

  return instance_idx + (function_idx * suite->instances.size()) +
         (dimension_idx * suite->instances.size() * suite->functions.size());
}

/**
 * @param suite The suite.
 * @param problem_index Index of the problem in the suite (starting with 0).
 * @param function_idx Pointer to the index of the function, which is set by
 * this function.
 * @param dimension_idx Pointer to the index of the dimension, which is set by
 * this function.
 * @param instance_idx Pointer to the index of the instance, which is set by
 * this function.
 */
void coco_suite_decode_problem_index(const coco_suite_t *suite, const size_t problem_index, size_t *function_idx,
                                     size_t *dimension_idx, size_t *instance_idx) {

  if (problem_index > (suite->instances.size() * suite->functions.size() *
                       suite->dimensions.size()) -
                          1) {
    coco_warning("coco_suite_decode_problem_index(): problem_index too large");
    *function_idx = 0;
    *instance_idx = 0;
    *dimension_idx = 0;
    return;
  }

  *instance_idx = problem_index % suite->instances.size();
  *function_idx = (problem_index / suite->instances.size()) % suite->functions.size();
  *dimension_idx = problem_index / (suite->instances.size() * suite->functions.size());
}
