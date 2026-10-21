/**
 * Tests of the public C interface (coco.h) that use nothing else, compiled as C.
 *
 * The language bindings also use the functions declared below, which are not part of coco.h.
 */
#include <stdio.h>
#include <string.h>

#include "coco.h"
#include "minunit.h"

coco_problem_t *coco_get_bbob_problem(const size_t function, const size_t dimension, const size_t instance);
double coco_problem_get_best_value(const coco_problem_t *problem);
double coco_problem_get_final_target_fvalue1(const coco_problem_t *problem);
void coco_suite_forget_current_problem(coco_suite_t *suite);

/** The observers write into this folder (relative to the working directory), which is removed afterwards. */
static const char *outer_folder = "test_coco_api_exdata";

/**
 * Silences the info messages of the tests and removes the output of the observers.
 */
static void setup(void) {
  coco_set_log_level("warning");
  coco_remove_directory(outer_folder);
}

static void teardown(void) {
  coco_remove_directory(outer_folder);
  coco_set_log_level("info");
}

static int starts_with(const char *string, const char *prefix) { return strncmp(string, prefix, strlen(prefix)) == 0; }

MU_TEST(test_version) {
  mu_check(coco_version != NULL);
  mu_check(strlen(coco_version) > 0);
}

/**
 * Tests coco_set_log_level, which returns the previous level, and the message functions at all levels.
 */
MU_TEST(test_log_level_and_messages) {
  mu_check(strcmp(coco_set_log_level("debug"), "warning") == 0);
  mu_check(strcmp(coco_set_log_level("info"), "debug") == 0);
  mu_check(strcmp(coco_set_log_level("error"), "info") == 0);
  /* An empty level only returns the current one */
  mu_check(strcmp(coco_set_log_level(""), "error") == 0);
  mu_check(strcmp(coco_set_log_level(NULL), "error") == 0);
  /* An unknown level leaves the level unchanged */
  coco_set_log_level("no such level");
  mu_check(strcmp(coco_set_log_level("warning"), "error") == 0);

  /* Below the log level, nothing is printed */
  coco_debug("test_coco_api: debug message %d", 1);
  coco_info("test_coco_api: info message %d", 2);
  coco_info_partial("test_coco_api: partial info message %d", 3);
  /* At the log level, the messages are printed */
  coco_set_log_level("debug");
  coco_debug("test_coco_api: debug message %s", "(expected)");
  coco_info("test_coco_api: info message %s", "(expected)");
  coco_info_partial("test_coco_api: partial info message %s\n", "(expected)");
  coco_warning("test_coco_api: warning %s", "(expected)");
  coco_set_log_level("warning");
}

MU_TEST(test_memory) {
  double *vector = coco_allocate_vector(5);
  char *memory = (char *)coco_allocate_memory(16);
  size_t i;

  for (i = 0; i < 5; ++i)
    vector[i] = (double)i;
  mu_check(vector[4] == 4.0);
  memset(memory, 'x', 16);
  mu_check(memory[15] == 'x');
  coco_free_memory(vector);
  coco_free_memory(memory);
}

MU_TEST(test_strdupf) {
  char *string = coco_strdupf("f%02d_%s_%.1f", 7, "bbob", 0.5);
  mu_check(strcmp(string, "f07_bbob_0.5") == 0);
  coco_free_memory(string);

  string = coco_strdupf("%s", "");
  mu_check(strcmp(string, "") == 0);
  coco_free_memory(string);
}

MU_TEST(test_random) {
  coco_random_state_t *random1 = coco_random_new(4711);
  coco_random_state_t *random2 = coco_random_new(4711);
  coco_random_state_t *random3 = coco_random_new(4712);
  double sum = 0, sum_of_squares = 0, value;
  int same = 1, differs = 0, in_range = 1;
  const int n = 20000;
  int i;

  /* The same seed gives the same numbers, beyond the size of the generator's buffer */
  for (i = 0; i < 1000; ++i) {
    value = coco_random_uniform(random1);
    same = same && (value == coco_random_uniform(random2));
    differs = differs || (value != coco_random_uniform(random3));
    in_range = in_range && (value >= 0.0) && (value < 1.0);
  }
  mu_check(same);
  mu_check(differs);
  mu_check(in_range);

  for (i = 0; i < n; ++i) {
    value = coco_random_normal(random1);
    sum += value;
    sum_of_squares += value * value;
  }
  mu_check(fabs(sum / n) < 0.05);
  mu_check(fabs(sum_of_squares / n - 1.0) < 0.05);

  coco_random_free(random1);
  coco_random_free(random2);
  coco_random_free(random3);
}

MU_TEST(test_suite_size_and_lookup) {
  coco_suite_t *suite = coco_suite("bbob", "instances: 1-3", "dimensions: 2,5");
  coco_problem_t *problem;

  mu_check(coco_suite_get_number_of_problems(suite) == 24 * 6 * 3);

  problem = coco_suite_get_problem_by_function_dimension_instance(suite, 3, 5, 2);
  mu_check(problem != NULL);
  mu_check(strcmp(coco_problem_get_id(problem), "bbob_f003_i02_d05") == 0);
  mu_check(coco_problem_get_suite_dep_index(problem) == coco_suite_encode_problem_index(suite, 2, 2, 1));
  coco_problem_free(problem);

  /* Function, dimension or instance not in the suite */
  mu_check(coco_suite_get_problem_by_function_dimension_instance(suite, 25, 5, 2) == NULL);
  mu_check(coco_suite_get_problem_by_function_dimension_instance(suite, 3, 4, 2) == NULL);
  mu_check(coco_suite_get_problem_by_function_dimension_instance(suite, 3, 5, 4) == NULL);
  coco_suite_free(suite);
}

/**
 * Tests that coco_suite_get_next_problem iterates by instances, then functions, then dimensions.
 */
MU_TEST(test_suite_iteration) {
  coco_suite_t *suite = coco_suite("bbob", "instances: 1-2", "dimensions: 2,3 function_indices: 1,2");
  coco_problem_t *problem;
  const char *expected[] = {"bbob_f001_i01_d02", "bbob_f001_i02_d02", "bbob_f002_i01_d02", "bbob_f002_i02_d02",
                            "bbob_f001_i01_d03", "bbob_f001_i02_d03", "bbob_f002_i01_d03", "bbob_f002_i02_d03"};
  size_t count = 0;

  while ((problem = coco_suite_get_next_problem(suite, NULL)) != NULL) {
    mu_check(count < 8);
    if (count >= 8)
      break;
    mu_check(strcmp(coco_problem_get_id(problem), expected[count]) == 0);
    count++;
  }
  mu_check(count == 8);
  coco_suite_free(suite);
}

/**
 * Tests that the caller owns the current problem after coco_suite_forget_current_problem.
 */
MU_TEST(test_suite_forget_current_problem) {
  coco_suite_t *suite = coco_suite("bbob", "instances: 1", "dimensions: 2 function_indices: 1");
  coco_problem_t *problem = coco_suite_get_next_problem(suite, NULL);
  double x[2] = {0, 0}, y;

  coco_suite_forget_current_problem(suite);
  coco_suite_free(suite);
  coco_evaluate_function(problem, x, &y);
  mu_check(coco_problem_get_evaluations(problem) == 1);
  coco_problem_free(problem);
}

MU_TEST(test_problem_properties) {
  coco_problem_t *problem = coco_get_bbob_problem(1, 2, 1);
  const double *lower = coco_problem_get_smallest_values_of_interest(problem);
  const double *upper = coco_problem_get_largest_values_of_interest(problem);
  double initial[2];

  mu_check(strcmp(coco_problem_get_id(problem), "bbob_f001_i01_d02") == 0);
  mu_check(strlen(coco_problem_get_name(problem)) > 0);
  mu_check(strlen(coco_problem_get_type(problem)) > 0);
  mu_check(coco_problem_get_dimension(problem) == 2);
  mu_check(coco_problem_get_number_of_objectives(problem) == 1);
  mu_check(coco_problem_get_number_of_constraints(problem) == 0);
  mu_check(coco_problem_get_number_of_integer_variables(problem) == 0);
  mu_check(lower[0] == -5.0 && lower[1] == -5.0);
  mu_check(upper[0] == 5.0 && upper[1] == 5.0);
  /* Without a given initial solution, it is the center of the region of interest */
  coco_problem_get_initial_solution(problem, initial);
  mu_check(initial[0] == 0.0 && initial[1] == 0.0);
  mu_check(coco_problem_get_final_target_fvalue1(problem) == coco_problem_get_best_value(problem) + 1e-8);
  coco_problem_free(problem);
}

/**
 * Tests the evaluation counters, the best observed value and the handling of NaN and infinite input.
 */
MU_TEST(test_problem_evaluation) {
  coco_problem_t *problem = coco_get_bbob_problem(1, 2, 1);
  const double best_value = coco_problem_get_best_value(problem);
  double x[2] = {1.0, 1.0}, y, y_far;

  mu_check(coco_problem_get_evaluations(problem) == 0);
  mu_check(coco_problem_final_target_hit(problem) == 0);

  coco_evaluate_function(problem, x, &y);
  mu_check(y > best_value);
  mu_check(coco_problem_get_evaluations(problem) == 1);
  mu_check(coco_problem_get_best_observed_fvalue1(problem) == y);

  /* A worse point does not change the best observed value */
  x[0] = 4.0;
  coco_evaluate_function(problem, x, &y_far);
  mu_check(y_far > y);
  mu_check(coco_problem_get_evaluations(problem) == 2);
  mu_check(coco_problem_get_best_observed_fvalue1(problem) == y);

  /* NaN and infinite decision values are not evaluated and not counted */
  x[0] = NAN;
  coco_evaluate_function(problem, x, &y);
  mu_check(isnan(y));
  x[0] = -INFINITY;
  coco_evaluate_function(problem, x, &y);
  mu_check(y == INFINITY);
  mu_check(coco_problem_get_evaluations(problem) == 2);
  mu_check(coco_problem_get_evaluations_constraints(problem) == 0);
  mu_check(coco_problem_final_target_hit(problem) == 0);

  /* An unobserved problem does not implement recommendations (only a warning) */
  coco_set_log_level("error");
  coco_recommend_solution(problem, x);
  coco_set_log_level("warning");
  coco_problem_free(problem);
}

/**
 * Tests that the best parameter is the optimum and that retrieving it taints the problem.
 */
MU_TEST(test_problem_best_parameter) {
  coco_problem_t *problem = coco_get_bbob_problem(1, 3, 2);
  double best_parameter[3], y;

  mu_check(coco_problem_is_tainted(problem) == 0);
  mu_check(coco_problem_get_best_parameter(problem, best_parameter) == 1);
  mu_check(coco_problem_is_tainted(problem) == 1);

  coco_evaluate_function(problem, best_parameter, &y);
  mu_check(y == coco_problem_get_best_value(problem));
  mu_check(coco_problem_get_best_observed_fvalue1(problem) == y);
  mu_check(coco_problem_final_target_hit(problem) == 1);

  /* Writes the best parameter to a file in the working directory */
  bbob_problem_best_parameter_print(problem);
  mu_check(remove("._bbob_problem_best_parameter.txt") == 0);
  coco_problem_free(problem);
}

MU_TEST(test_biobj_problem) {
  coco_suite_t *suite = coco_suite("bbob-biobj", "instances: 1", "dimensions: 2 function_indices: 1");
  coco_problem_t *problem = coco_suite_get_problem(suite, 0);
  const double *nadir = coco_problem_get_largest_fvalues_of_interest(problem);
  double x[2] = {0, 0}, y[2], best_parameter[2];
  FILE *file;

  mu_check(coco_problem_get_number_of_objectives(problem) == 2);
  mu_check(isfinite(nadir[0]) && isfinite(nadir[1]));
  coco_evaluate_function(problem, x, y);
  mu_check(isfinite(y[0]) && isfinite(y[1]));
  /* Only single-objective problems have a best parameter */
  mu_check(coco_problem_get_best_parameter(problem, best_parameter) == 0);
  mu_check(coco_problem_is_tainted(problem) == 0);

  /* Writes the best parameters of both objectives to a file in the working directory */
  bbob_biobj_problem_best_parameter_print(problem);
  file = fopen("._bbob_biobj_problem_best_parameter.txt", "r");
  mu_check(file != NULL);
  if (file != NULL) {
    double value;
    int count = 0;
    while (fscanf(file, "%lf", &value) == 1)
      count++;
    mu_check(count == 4);
    fclose(file);
  }
  mu_check(remove("._bbob_biobj_problem_best_parameter.txt") == 0);
  coco_problem_free(problem);
  coco_suite_free(suite);
}

/**
 * Tests the constraint evaluation and that the initial solution of a constrained problem is feasible.
 */
MU_TEST(test_constrained_problem) {
  coco_suite_t *suite = coco_suite("bbob-constrained", "instances: 1", "dimensions: 2 function_indices: 1");
  coco_problem_t *problem = coco_suite_get_problem(suite, 0);
  const size_t number_of_constraints = coco_problem_get_number_of_constraints(problem);
  double x[2], y, *constraints = coco_allocate_vector(number_of_constraints);
  size_t i;

  mu_check(number_of_constraints > 0);
  coco_problem_get_initial_solution(problem, x);
  coco_evaluate_constraint(problem, x, constraints);
  for (i = 0; i < number_of_constraints; ++i)
    mu_check(constraints[i] <= 0.0);
  mu_check(coco_problem_get_evaluations_constraints(problem) == 1);
  coco_evaluate_function(problem, x, &y);
  mu_check(coco_problem_get_evaluations(problem) == 1);
  mu_check(coco_problem_get_evaluations_constraints(problem) == 1);

  /* NaN decision values give NaN constraints and are not counted */
  x[1] = NAN;
  coco_evaluate_constraint(problem, x, constraints);
  mu_check(isnan(constraints[0]));
  mu_check(coco_problem_get_evaluations_constraints(problem) == 1);

  coco_free_memory(constraints);
  coco_problem_free(problem);
  coco_suite_free(suite);
}

/**
 * Tests that the integer variables come first and that the initial solution is integral in them.
 */
MU_TEST(test_mixint_problem) {
  coco_suite_t *suite = coco_suite("bbob-mixint", "instances: 1", "dimensions: 5 function_indices: 1");
  coco_problem_t *problem = coco_suite_get_problem(suite, 0);
  const size_t number_of_integer_variables = coco_problem_get_number_of_integer_variables(problem);
  const double *lower = coco_problem_get_smallest_values_of_interest(problem);
  const double *upper = coco_problem_get_largest_values_of_interest(problem);
  double x[5];
  size_t i;

  mu_check(number_of_integer_variables > 0);
  mu_check(number_of_integer_variables < 5);
  coco_problem_get_initial_solution(problem, x);
  for (i = 0; i < 5; ++i)
    mu_check(lower[i] <= x[i] && x[i] <= upper[i]);
  for (i = 0; i < number_of_integer_variables; ++i)
    mu_check(x[i] == floor(x[i]));
  coco_problem_free(problem);
  coco_suite_free(suite);
}

/**
 * Tests that coco_reset_seeds makes the noise of the noisy suite reproducible.
 */
MU_TEST(test_reset_seeds) {
  coco_suite_t *suite = coco_suite("bbob-noisy", "instances: 1", "dimensions: 2 function_indices: 1");
  coco_problem_t *problem = coco_suite_get_problem(suite, 0);
  double x[2] = {1.0, 2.0}, y1, y2, y3;

  coco_reset_seeds();
  coco_evaluate_function(problem, x, &y1);
  coco_evaluate_function(problem, x, &y2);
  coco_reset_seeds();
  coco_evaluate_function(problem, x, &y3);
  mu_check(y1 != y2);
  mu_check(y1 == y3);
  coco_problem_free(problem);
  coco_suite_free(suite);
}

/**
 * Tests adding and removing an observer and the observer functions on an observed problem.
 */
MU_TEST(test_observer) {
  char *options = coco_strdupf("outer_folder: %s result_folder: test_observer", outer_folder);
  coco_observer_t *observer = coco_observer("bbob", options);
  /* The bbob logger needs the suite of the problem */
  coco_suite_t *suite = coco_suite("bbob", "instances: 1", "dimensions: 3");
  coco_problem_t *problem = coco_suite_get_problem_by_function_dimension_instance(suite, 2, 3, 1);
  coco_problem_t *observed;
  const char *name = coco_problem_get_name(problem);
  char *inner_name = coco_strdupf("%s", name);
  double x[3] = {1, 2, 3}, y;

  mu_check(observer != NULL);
  mu_check(strstr(coco_observer_get_result_folder(observer), "test_observer") != NULL);

  observed = coco_problem_add_observer(problem, observer);
  mu_check(starts_with(coco_problem_get_name(observed), "bbob("));
  coco_evaluate_function(observed, x, &y);
  coco_recommend_solution(observed, x);
  coco_observer_signal_restart(observer, observed);
  coco_evaluate_function(observed, x, &y);
  mu_check(coco_problem_get_evaluations(observed) == 2);

  problem = coco_problem_remove_observer(observed, observer);
  mu_check(strcmp(coco_problem_get_name(problem), inner_name) == 0);
  coco_problem_free(problem);
  coco_observer_free(observer);
  coco_suite_free(suite);
  coco_free_memory(inner_name);
  coco_free_memory(options);
}

/**
 * Tests the observer functions without an observer or with an inactive one; they only warn.
 */
MU_TEST(test_no_observer) {
  coco_problem_t *problem = coco_get_bbob_problem(1, 2, 1);
  char *options = coco_strdupf("outer_folder: %s log_nondominated: none compute_indicators: 0", outer_folder);
  coco_observer_t *inactive = coco_observer("bbob-biobj", options);

  mu_check(inactive != NULL);
  mu_check(coco_observer("no_observer", "") == NULL);

  coco_set_log_level("error");
  mu_check(coco_observer(NULL, "") == NULL);
  mu_check(coco_observer("", "") == NULL);
  mu_check(coco_observer("no_such_observer", options) == NULL);
  mu_check(strcmp(coco_observer_get_result_folder(NULL), "") == 0);
  mu_check(strcmp(coco_observer_get_result_folder(inactive), "") == 0);
  mu_check(coco_problem_add_observer(problem, NULL) == problem);
  mu_check(coco_problem_add_observer(problem, inactive) == problem);
  mu_check(coco_problem_remove_observer(problem, NULL) == problem);
  mu_check(coco_problem_remove_observer(problem, inactive) == problem);
  coco_observer_signal_restart(NULL, problem);
  coco_observer_signal_restart(inactive, problem);
  coco_set_log_level("warning");
  mu_check(coco_problem_add_observer(NULL, inactive) == NULL);

  coco_observer_free(inactive);
  coco_problem_free(problem);
  coco_free_memory(options);
}

/**
 * Tests that the bi-objective logger in "read" mode keeps only the nondominated fed solutions.
 */
MU_TEST(test_logger_biobj_feed_solution) {
  char *options = coco_strdupf("outer_folder: %s result_folder: test_feed log_nondominated: read", outer_folder);
  coco_observer_t *observer = coco_observer("bbob-biobj", options);
  coco_suite_t *suite = coco_suite("bbob-biobj", "instances: 1", "dimensions: 2 function_indices: 1");
  coco_problem_t *problem = coco_suite_get_next_problem(suite, observer);
  const double *nadir = coco_problem_get_largest_fvalues_of_interest(problem);
  double y[2];

  y[0] = nadir[0] - 1.0;
  y[1] = nadir[1] - 1.0;
  mu_check(coco_logger_biobj_feed_solution(problem, 1, y) == 1);
  /* Dominated by the first solution */
  y[0] = nadir[0] - 0.5;
  y[1] = nadir[1] - 0.5;
  mu_check(coco_logger_biobj_feed_solution(problem, 2, y) == 0);
  /* Not dominated */
  y[0] = nadir[0] - 2.0;
  y[1] = nadir[1] - 0.75;
  mu_check(coco_logger_biobj_feed_solution(problem, 3, y) == 1);

  coco_suite_free(suite);
  coco_observer_free(observer);
  coco_free_memory(options);
}

MU_TEST(test_archive) {
  coco_archive_t *archive = coco_archive("bbob-biobj", 1, 2, 1);
  const size_t extremes = coco_archive_get_number_of_solutions(archive);
  const char *text;
  size_t count = 0;

  /* The archive starts with the two extreme solutions */
  mu_check(extremes == 2);
  mu_check(coco_archive_get_hypervolume(archive) == 0.0);
  while (strlen(text = coco_archive_get_next_solution_text(archive)) > 0)
    count++;
  mu_check(count == 2);
  coco_archive_free(archive);
}

MU_TEST(test_remove_directory) {
  /* A directory that does not exist is not an error */
  mu_check(coco_remove_directory("test_coco_api_no_such_directory") == 0);
}

MU_TEST_SUITE(test_all_coco_api) {
  MU_SUITE_CONFIGURE(setup, teardown);
  MU_RUN_TEST(test_version);
  MU_RUN_TEST(test_log_level_and_messages);
  MU_RUN_TEST(test_memory);
  MU_RUN_TEST(test_strdupf);
  MU_RUN_TEST(test_random);
  MU_RUN_TEST(test_suite_size_and_lookup);
  MU_RUN_TEST(test_suite_iteration);
  MU_RUN_TEST(test_suite_forget_current_problem);
  MU_RUN_TEST(test_problem_properties);
  MU_RUN_TEST(test_problem_evaluation);
  MU_RUN_TEST(test_problem_best_parameter);
  MU_RUN_TEST(test_biobj_problem);
  MU_RUN_TEST(test_constrained_problem);
  MU_RUN_TEST(test_mixint_problem);
  MU_RUN_TEST(test_reset_seeds);
  MU_RUN_TEST(test_observer);
  MU_RUN_TEST(test_no_observer);
  MU_RUN_TEST(test_logger_biobj_feed_solution);
  MU_RUN_TEST(test_archive);
  MU_RUN_TEST(test_remove_directory);
}

int main(void) {
  MU_RUN_SUITE(test_all_coco_api);
  MU_REPORT();
  return MU_EXIT_CODE;
}
