/**
 * @file suite_biobj_utilities.cpp
 * @brief Implementation of some functions (mostly handling instances) used by the bi-objective suites.
 *
 * @note Because some bi-objective problems constructed from two single-objective ones have a single optimal
 * value, some care must be taken when selecting the instances. The already verified instances are stored in
 * suite_biobj_instances. If a new instance of the problem is called, a check ensures that the two underlying
 * single-objective instances create a true bi-objective problem. However, these new instances need to be
 * manually added to suite_biobj_instances, otherwise they will be computed each time the suite constructor
 * is invoked with these instances.
 */

#include <array>

#include "coco.h"
#include "suite_biobj_best_values_hyp.cpp"

/**
 * @brief The array of triples biobj_instance - problem1_instance - problem2_instance connecting bi-objective
 * suite instances to the instances of the bbob suite.
 *
 * It should be updated with new instances when/if they are chosen.
 */
static const size_t suite_biobj_instances[][3] = {{1, 2, 4},    {2, 3, 5},    {3, 7, 8},    {4, 9, 10},   {5, 11, 12},
                                                  {6, 13, 14},  {7, 15, 16},  {8, 17, 18},  {9, 19, 21},  {10, 21, 22},
                                                  {11, 23, 24}, {12, 25, 26}, {13, 27, 28}, {14, 29, 30}, {15, 31, 34}};

/**
 * @brief The new instances (triples equal in form to suite_biobj_instances) that are computed when an
 * instance that is not in suite_biobj_instances is being invoked.
 */
using suite_biobj_new_instances_t = std::vector<std::array<size_t, 3>>;

/**
 * @brief  Performs a few checks and returns whether the two given problem instances should break the search
 * for new instances in suite_biobj_get_new_instance().
 */
static bool suite_biobj_check_inst_consistency(const size_t dimension, size_t function1, size_t instance1,
                                               size_t function2, size_t instance2) {
  bool break_search = false;
  const double apart_enough = 1e-4;

  coco_problem_ptr_t problem1 = coco_make_bbob_problem(function1, dimension, instance1);
  coco_problem_ptr_t problem2 = coco_make_bbob_problem(function2, dimension, instance2);

  /* Set smallest and largest values of interest to some value (not important which, it just needs to be a
   * vector of doubles of the right dimension) */
  const std::vector<double> smallest_values_of_interest(dimension, -100);
  const std::vector<double> largest_values_of_interest(dimension, 100);
  const coco_problem_t *problem1_view = problem1.get();
  const coco_problem_t *problem2_view = problem2.get();
  const coco_problem_ptr_t problem = std::make_unique<coco_problem_stacked_t>(
      std::move(problem1), std::move(problem2), smallest_values_of_interest, largest_values_of_interest);

  /* Check whether the ideal and nadir points are too close in the objective space */
  double norm = mo_get_norm(problem->get_best_value(), problem->get_nadir_value());
  if (norm < 1e-1) { /* TODO How to set this value in a sensible manner? */
    coco_debug(
        "suite_biobj_check_inst_consistency(): The ideal and nadir points of %s are too close in the objective space",
        problem->get_id().c_str());
    coco_debug("norm = %e, ideal = %e\t%e, nadir = %e\t%e", norm, problem->get_best_value()[0],
               problem->get_best_value()[1], problem->get_nadir_value()[0], problem->get_nadir_value()[1]);
    break_search = true;
  }

  /* Check whether the extreme optimal points are too close in the decision space */
  norm = mo_get_norm(problem1_view->get_best_parameter(), problem2_view->get_best_parameter());
  if (norm < apart_enough) {
    coco_debug("suite_biobj_check_inst_consistency(): The extreme points of %s are too close in the decision space",
               problem->get_id().c_str());
    coco_debug("norm = %e", norm);
    break_search = true;
  }

  return break_search;
}

/**
 * @brief Computes the instance number of the second problem/objective so that the resulting bi-objective
 * problem has more than a single optimal solution.
 *
 * Starts by setting instance2 = instance1 + 1 and increases this number until an appropriate instance has
 * been found (or until a maximum number of tries has been reached, in which case it throws a coco_error).
 * An appropriate instance is the one for which the resulting bi-objective problem (in any considered
 * dimension) has the ideal and nadir points apart enough in the objective space and the extreme optimal
 * points apart enough in the decision space. When the instance has been found, it is output through
 * coco_warning, so that the user can see it and eventually manually add it to suite_biobj_instances.
 */
static size_t suite_biobj_get_new_instance(suite_biobj_new_instances_t &new_instances, const size_t instance,
                                           const size_t instance1, const size_t *bbob_functions,
                                           const size_t num_bbob_functions, const size_t *sel_bbob_functions,
                                           const size_t num_sel_bbob_functions, const size_t *dimensions,
                                           const size_t num_dimensions) {

  size_t instance2 = 0;
  size_t num_tries = 0;
  const size_t max_tries = 1000;
  bool appropriate_instance_found = false, break_search, warning_produced = false;
  size_t function1, function2, dimension;

  while ((!appropriate_instance_found) && (num_tries < max_tries)) {
    num_tries++;
    instance2 = instance1 + num_tries;
    break_search = false;

    /* An instance is "appropriate" if the ideal and nadir points in the objective space and the two
     * extreme optimal points in the decisions space are apart enough for all problems (all dimensions
     * and function combinations); therefore iterate over all dimensions and function combinations */

    for (size_t f1 = 0; (f1 < num_bbob_functions - 1) && !break_search; f1++) {
      function1 = bbob_functions[f1];
      for (size_t f2 = f1 + 1; (f2 < num_bbob_functions) && !break_search; f2++) {
        function2 = bbob_functions[f2];
        for (size_t d = 0; (d < num_dimensions) && !break_search; d++) {
          dimension = dimensions[d];

          if (dimension == 0) {
            if (!warning_produced)
              coco_warning("suite_biobj_get_new_instance(): remove filtering of dimensions to get generally acceptable "
                           "instances!");
            warning_produced = true;
            continue;
          }

          break_search = suite_biobj_check_inst_consistency(dimension, function1, instance1, function2, instance2);
        }
      }
    }

    /* Finally, check all functions (f,f) with f in {f1, f2, f6, f8, f13, f14, f15, f17, f20, f21}: */
    for (size_t f1 = 0; (f1 < num_sel_bbob_functions) && !break_search; f1++) {
      function1 = sel_bbob_functions[f1];
      function2 = sel_bbob_functions[f1];
      for (size_t d = 0; (d < num_dimensions) && !break_search; d++) {
        dimension = dimensions[d];

        if (dimension == 0) {
          if (!warning_produced)
            coco_warning("suite_biobj_get_new_instance(): remove filtering of dimensions to get generally acceptable "
                         "instances!");
          warning_produced = true;
          continue;
        }

        break_search = suite_biobj_check_inst_consistency(dimension, function1, instance1, function2, instance2);
      }
    }

    if (break_search) {
      /* The search was broken, continue with next instance2 */
      continue;
    } else {
      /* An appropriate instance was found */
      appropriate_instance_found = true;
      coco_info("suite_biobj_get_new_instance(): Instance %zu created from instances %zu and %zu", instance, instance1,
                instance2);

      /* Save the instance to new_instances */
      new_instances.push_back({instance, instance1, instance2});
    }
  }

  if (!appropriate_instance_found) {
    coco_error("suite_biobj_get_new_instance(): Could not find suitable instance %zu in %zu tries", instance,
               num_tries);
    return 0; /* Never reached */
  }

  return instance2;
}

/**
 * @brief Creates and returns a bi-objective problem without needing a suite.
 *
 * Useful for creating suites based on the bi-objective problems.
 *
 * Creates the bi-objective problem by constructing it from two single-objective problems. If the
 * invoked instance number is not in suite_biobj_instances, the function uses the following formula
 * to construct a new appropriate instance:
 *   problem1_instance = 2 * biobj_instance + 1
 *   problem2_instance = problem1_instance + 1
 *
 * If needed, problem2_instance is increased (see also the explanation in suite_biobj_get_new_instance).
 *
 * @param function Function
 * @param dimension Dimension
 * @param instance Instance
 * @param coco_get_problem_function The function that is used to access the single-objective problem.
 * @param new_instances The new instances computed so far (new ones are added).
 * @param dimensions An array of dimensions to take into account when creating new instances.
 * @param num_dimensions The number of dimensions to take into account when creating new instances.
 * @return The problem that corresponds to the given parameters.
 */
static coco_problem_ptr_t coco_get_biobj_problem(const size_t function, const size_t dimension, const size_t instance,
                                                 const coco_get_problem_function_t coco_get_problem_function,
                                                 suite_biobj_new_instances_t &new_instances, const size_t *dimensions,
                                                 const size_t num_dimensions) {

  /* Selected functions from the bbob suite that are used to construct the original bbob-biobj suite. */
  const size_t sel_bbob_functions[] = {1, 2, 6, 8, 13, 14, 15, 17, 20, 21};
  const size_t num_sel_bbob_functions = 10;
  /* All functions from the bbob suite for later use during instance generation. */
  const size_t all_bbob_functions[] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                                       13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
  const size_t num_all_bbob_functions = 24;

  coco_problem_ptr_t problem1;
  coco_problem_ptr_t problem2;
  coco_problem_ptr_t problem;
  size_t instance1 = 0, instance2 = 0;
  size_t function1_idx, function2_idx;
  const size_t function_idx = function - 1;

  const size_t num_existing_instances = sizeof(suite_biobj_instances) / sizeof(suite_biobj_instances[0]);
  bool instance_found = false;

  const std::vector<double> smallest_values_of_interest(dimension, -100);
  const std::vector<double> largest_values_of_interest(dimension, 100);

  /* Determine the corresponding single-objective function indices */
  if (function_idx < 55) {
    /* A "magic" formula to compute the BBOB function index from the bi-objective function index */
    function1_idx = num_sel_bbob_functions -
                    coco_double_to_size_t(floor(-0.5 + sqrt(0.25 + 2.0 * static_cast<double>(55 - function_idx - 1)))) -
                    1;
    function2_idx = function_idx - (function1_idx * num_sel_bbob_functions) + (function1_idx * (function1_idx + 1)) / 2;
  } else if (function_idx == 55) { /* There is not a simple "magic" formula for functions >= 55 */
    function1_idx = 0;
    function2_idx = 2;
  } else if (function_idx == 56) {
    function1_idx = 0;
    function2_idx = 3;
  } else if (function_idx == 57) {
    function1_idx = 0;
    function2_idx = 4;
  } else if (function_idx == 58) {
    function1_idx = 1;
    function2_idx = 2;
  } else if (function_idx == 59) {
    function1_idx = 1;
    function2_idx = 3;
  } else if (function_idx == 60) {
    function1_idx = 1;
    function2_idx = 4;
  } else if (function_idx == 61) {
    function1_idx = 2;
    function2_idx = 3;
  } else if (function_idx == 62) {
    function1_idx = 2;
    function2_idx = 4;
  } else if (function_idx == 63) {
    function1_idx = 3;
    function2_idx = 4;
  } else if (function_idx == 64) {
    function1_idx = 5;
    function2_idx = 6;
  } else if (function_idx == 65) {
    function1_idx = 5;
    function2_idx = 8;
  } else if (function_idx == 66) {
    function1_idx = 6;
    function2_idx = 7;
  } else if (function_idx == 67) {
    function1_idx = 6;
    function2_idx = 8;
  } else if (function_idx == 68) {
    function1_idx = 7;
    function2_idx = 8;
  } else if (function_idx == 69) {
    function1_idx = 9;
    function2_idx = 10;
  } else if (function_idx == 70) {
    function1_idx = 9;
    function2_idx = 11;
  } else if (function_idx == 71) {
    function1_idx = 9;
    function2_idx = 12;
  } else if (function_idx == 72) {
    function1_idx = 9;
    function2_idx = 13;
  } else if (function_idx == 73) {
    function1_idx = 10;
    function2_idx = 11;
  } else if (function_idx == 74) {
    function1_idx = 10;
    function2_idx = 12;
  } else if (function_idx == 75) {
    function1_idx = 10;
    function2_idx = 13;
  } else if (function_idx == 76) {
    function1_idx = 11;
    function2_idx = 12;
  } else if (function_idx == 77) {
    function1_idx = 11;
    function2_idx = 13;
  } else if (function_idx == 78) {
    function1_idx = 14;
    function2_idx = 17;
  } else if (function_idx == 79) {
    function1_idx = 14;
    function2_idx = 18;
  } else if (function_idx == 80) {
    function1_idx = 16;
    function2_idx = 17;
  } else if (function_idx == 81) {
    function1_idx = 16;
    function2_idx = 18;
  } else if (function_idx == 82) {
    function1_idx = 17;
    function2_idx = 18;
  } else if (function_idx == 83) {
    function1_idx = 19;
    function2_idx = 21;
  } else if (function_idx == 84) {
    function1_idx = 19;
    function2_idx = 22;
  } else if (function_idx == 85) {
    function1_idx = 19;
    function2_idx = 23;
  } else if (function_idx == 86) {
    function1_idx = 20;
    function2_idx = 21;
  } else if (function_idx == 87) {
    function1_idx = 20;
    function2_idx = 22;
  } else if (function_idx == 88) {
    function1_idx = 20;
    function2_idx = 23;
  } else if (function_idx == 89) {
    function1_idx = 21;
    function2_idx = 22;
  } else if (function_idx == 90) {
    function1_idx = 21;
    function2_idx = 23;
  } else if (function_idx == 91) {
    function1_idx = 22;
    function2_idx = 23;
  } else {
    coco_error("coco_get_biobj_problem(): Invalid function index %zu.", function_idx);
  }

  /* Determine the instances */

  /* First search for the instance in suite_biobj_instances */
  for (size_t i = 0; i < num_existing_instances; i++) {
    if (suite_biobj_instances[i][0] == instance) {
      /* The instance has been found in suite_biobj_instances */
      instance1 = suite_biobj_instances[i][1];
      instance2 = suite_biobj_instances[i][2];
      instance_found = true;
      break;
    }
  }

  if (!instance_found) {
    /* Next, search for instance in new_instances */
    for (const auto &new_instance : new_instances) {
      if (new_instance[0] == instance) {
        /* The instance has been found in new_instances */
        instance1 = new_instance[1];
        instance2 = new_instance[2];
        instance_found = true;
        break;
      }
    }
  }

  if (!instance_found) {
    /* Finally, if the instance is not found, create a new one */
    /* A simple formula to set the first instance */
    instance1 = 2 * instance + 1;
    instance2 = suite_biobj_get_new_instance(new_instances, instance, instance1, all_bbob_functions,
                                             num_all_bbob_functions, sel_bbob_functions, num_sel_bbob_functions,
                                             dimensions, num_dimensions);
  }

  /* Construct the problem based on the function index and dimension */
  if (function_idx < 55) {
    problem1 = coco_get_problem_function(sel_bbob_functions[function1_idx], dimension, instance1);
    problem2 = coco_get_problem_function(sel_bbob_functions[function2_idx], dimension, instance2);
    /* Store function numbers of the underlying problems */
    problem1->set_suite_dep_function(sel_bbob_functions[function1_idx]);
    problem2->set_suite_dep_function(sel_bbob_functions[function2_idx]);

  } else {
    problem1 = coco_get_problem_function(all_bbob_functions[function1_idx], dimension, instance1);
    problem2 = coco_get_problem_function(all_bbob_functions[function2_idx], dimension, instance2);
    problem1->set_suite_dep_function(all_bbob_functions[function1_idx]);
    problem2->set_suite_dep_function(all_bbob_functions[function2_idx]);
  }

  const std::string problem_type1 = problem1->get_type();
  const std::string problem_type2 = problem2->get_type();
  problem = std::make_unique<coco_problem_stacked_t>(std::move(problem1), std::move(problem2),
                                                     smallest_values_of_interest, largest_values_of_interest);

  /* Use the standard stacked problem_id as problem_name and construct a new problem_id */
  problem->set_name("%s", problem->get_id().c_str());
  /* Attention! Any change to the problem id affects also archive processing! */
  problem->set_id("bbob-biobj_f%02zu_i%02zu_d%02zu", function, instance, dimension);

  /* Construct problem type */
  problem->set_type("%s_%s", problem_type1.c_str(), problem_type2.c_str());

  return problem;
}

/**
 * @brief Saves the best known value for the hypervolume indicator matching the given key.
 *
 * If a suite does not have known optima or it has known optima but the key is not found,
 * the default value is used.
 */
static void suite_biobj_get_best_hyp_value(const int known_optima, const std::string &key, double *value) {

  static const double default_value = 1.0;
  size_t count;
  *value = default_value;

  if (known_optima) {
    std::vector<char> curr_key;
    count = sizeof(suite_biobj_best_values_hyp) / sizeof(char *);
    for (size_t i = 0; i < count; i++) {
      curr_key.resize(strlen(suite_biobj_best_values_hyp[i]) + 1); /* The key cannot be longer than the entry */
      sscanf(suite_biobj_best_values_hyp[i], "%s %lf", curr_key.data(), value);
      if (key == curr_key.data()) {
        return;
      }
    }
    /* If it comes to this point, the key was not found */
    coco_warning("suite_biobj_get_best_hyp_value(): best value of %s could not be found; set to %f", key.c_str(),
                 default_value);
  } else {
    coco_warning("suite_biobj_get_best_hyp_value(): best value of %s is not known; set to %f", key.c_str(),
                 default_value);
  }
}
