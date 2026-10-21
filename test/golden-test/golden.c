/**
 * Golden master driver for COCO.
 *
 * Records the observable behaviour of COCO using only the public C interface so that refactorings of the
 * implementation can be checked for bit-exact equality.
 *
 * Modes:
 *
 *   golden values <suite> [--dump]
 *     Evaluates every problem of the suite (default instances) at a few deterministic points and prints one
 *     line per (function, dimension) pair containing a hash over all observable values. With --dump every
 *     value is printed instead (as hex float), which is useful to find out what changed.
 *
 *   golden observe <suite> <observer> <suite options>
 *     Runs a deterministic random search on the selected problems with the given observer attached so that
 *     the logger output can be compared.
 */

#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "coco.h"

/* Not part of coco.h but exported by the amalgamation and used by the Python bindings. */
double coco_problem_get_best_value(const coco_problem_t *problem);

#define NUMBER_OF_RANDOM_POINTS 3

/* --- Deterministic random numbers independent of COCO (splitmix64) ------------------------------------ */

static uint64_t rng_state;

static uint64_t rng_next(void) {
  uint64_t z = (rng_state += UINT64_C(0x9E3779B97F4A7C15));
  z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
  return z ^ (z >> 31);
}

/* Uniform in [0, 1) with 53 random bits */
static double rng_uniform(void) { return (double)(rng_next() >> 11) * 0x1.0p-53; }

/* --- Hashing (FNV-1a, 64 bit) --------------------------------------------------------------------------- */

static uint64_t hash;
static int dump;

static void hash_bytes(const void *data, size_t size) {
  const unsigned char *bytes = (const unsigned char *)data;
  size_t i;
  for (i = 0; i < size; ++i) {
    hash ^= bytes[i];
    hash *= UINT64_C(0x100000001B3);
  }
}

static void record_string(const char *label, const char *value) {
  if (value == NULL)
    value = "(null)";
  hash_bytes(value, strlen(value) + 1);
  if (dump)
    printf("  %s: %s\n", label, value);
}

static void record_size(const char *label, size_t value) {
  uint64_t v = (uint64_t)value;
  hash_bytes(&v, sizeof(v));
  if (dump)
    printf("  %s: %" PRIu64 "\n", label, v);
}

static void record_vector(const char *label, const double *values, size_t n) {
  size_t i;
  if (dump)
    printf("  %s:", label);
  for (i = 0; i < n; ++i) {
    double v = values[i];
    uint64_t bits;
    /* All NaNs are equal for our purposes, their payload and sign are not specified. */
    if (isnan(v))
      v = NAN;
    memcpy(&bits, &v, sizeof(bits));
    if (isnan(v))
      bits = UINT64_C(0x7FF8000000000000);
    hash_bytes(&bits, sizeof(bits));
    if (dump)
      printf(" %a", values[i]);
  }
  if (dump)
    printf("\n");
}

/* --- Values mode ---------------------------------------------------------------------------------------- */

static void evaluate_at(coco_problem_t *problem, const char *label, const double *x, double *y, double *c) {
  const size_t n = coco_problem_get_dimension(problem);
  const size_t m = coco_problem_get_number_of_objectives(problem);
  const size_t k = coco_problem_get_number_of_constraints(problem);

  if (dump)
    printf(" point %s\n", label);
  record_vector("x", x, n);
  coco_evaluate_function(problem, x, y);
  record_vector("f", y, m);
  if (k > 0) {
    coco_evaluate_constraint(problem, x, c);
    record_vector("c", c, k);
  }
}

static void record_problem(coco_problem_t *problem, int has_type) {
  const size_t n = coco_problem_get_dimension(problem);
  const size_t m = coco_problem_get_number_of_objectives(problem);
  const size_t k = coco_problem_get_number_of_constraints(problem);
  const double *lower = coco_problem_get_smallest_values_of_interest(problem);
  const double *upper = coco_problem_get_largest_values_of_interest(problem);
  double *x = (double *)malloc(n * sizeof(double));
  double *y = (double *)malloc(m * sizeof(double));
  double *c = (double *)malloc((k > 0 ? k : 1) * sizeof(double));
  double best_value;
  char label[32];
  size_t i, j;

  if (dump)
    printf("problem %s\n", coco_problem_get_id(problem));
  record_string("id", coco_problem_get_id(problem));
  record_string("name", coco_problem_get_name(problem));
  if (has_type) /* coco_problem_get_type() asserts that the type is set */
    record_string("type", coco_problem_get_type(problem));
  record_size("dimension", n);
  record_size("objectives", m);
  record_size("constraints", k);
  record_size("integer variables", coco_problem_get_number_of_integer_variables(problem));
  record_vector("lower", lower, n);
  record_vector("upper", upper, n);
  if (m == 1) {
    best_value = coco_problem_get_best_value(problem);
    record_vector("best value", &best_value, 1);
  } else {
    record_vector("largest f of interest", coco_problem_get_largest_fvalues_of_interest(problem), m);
  }

  coco_problem_get_initial_solution(problem, x);
  evaluate_at(problem, "initial", x, y, c);

  /* Seed per problem so that a change in one problem does not cascade into all following ones. */
  rng_state = coco_problem_get_suite_dep_index(problem);
  for (j = 0; j < NUMBER_OF_RANDOM_POINTS; ++j) {
    for (i = 0; i < n; ++i)
      x[i] = lower[i] + rng_uniform() * (upper[i] - lower[i]);
    sprintf(label, "random %lu", (unsigned long)j);
    evaluate_at(problem, label, x, y, c);
  }
  /* A point that is likely outside of the region of interest */
  for (i = 0; i < n; ++i)
    x[i] = lower[i] + (3.0 * rng_uniform() - 1.0) * (upper[i] - lower[i]);
  evaluate_at(problem, "outside", x, y, c);

  /* Points with a non-finite coordinate are handled before the problem is evaluated */
  coco_problem_get_initial_solution(problem, x);
  x[n / 2] = NAN;
  evaluate_at(problem, "nan", x, y, c);
  x[n / 2] = INFINITY;
  evaluate_at(problem, "+inf", x, y, c);
  x[n / 2] = -INFINITY;
  evaluate_at(problem, "-inf", x, y, c);

  /* Last because it taints the problem */
  if (coco_problem_get_best_parameter(problem, x))
    evaluate_at(problem, "best", x, y, c);

  record_size("evaluations", coco_problem_get_evaluations(problem));
  record_size("constraint evaluations", coco_problem_get_evaluations_constraints(problem));

  free(x);
  free(y);
  free(c);
}

static int values(const char *suite_name) {
  coco_suite_t *suite = coco_suite(suite_name, "", "");
  const size_t number_of_problems = coco_suite_get_number_of_problems(suite);
  size_t index, function_idx, dimension_idx, instance_idx;
  size_t group_function = 0, group_dimension = 0, group_problems = 0;

  coco_reset_seeds(); /* bbob-noisy uses a global random state */
  printf("suite %s problems %lu\n", suite_name, (unsigned long)number_of_problems);

  for (index = 0; index < number_of_problems; ++index) {
    coco_problem_t *problem = coco_suite_get_problem(suite, index);
    size_t function, dimension;

    coco_suite_decode_problem_index(suite, index, &function_idx, &dimension_idx, &instance_idx);
    function = coco_suite_get_function_from_function_index(suite, function_idx);
    dimension = coco_suite_get_dimension_from_dimension_index(suite, dimension_idx);

    if (group_problems > 0 && (function != group_function || dimension != group_dimension)) {
      if (!dump)
        printf("f%03lu d%03lu problems %lu hash %016" PRIx64 "\n", (unsigned long)group_function,
               (unsigned long)group_dimension, (unsigned long)group_problems, hash);
      group_problems = 0;
    }
    if (group_problems == 0) {
      hash = UINT64_C(0xCBF29CE484222325);
      group_function = function;
      group_dimension = dimension;
    }
    record_problem(problem, strcmp(suite_name, "toy") != 0);
    ++group_problems;
    coco_problem_free(problem);
  }
  if (group_problems > 0 && !dump)
    printf("f%03lu d%03lu problems %lu hash %016" PRIx64 "\n", (unsigned long)group_function,
           (unsigned long)group_dimension, (unsigned long)group_problems, hash);

  coco_suite_free(suite);
  return EXIT_SUCCESS;
}

/* --- Observe mode --------------------------------------------------------------------------------------- */

static int observe(const char *suite_name, const char *observer_name, const char *suite_options) {
  coco_suite_t *suite = coco_suite(suite_name, "", suite_options);
  coco_observer_t *observer = coco_observer(observer_name, "result_folder: golden algorithm_name: golden "
                                                           "algorithm_info: \"deterministic random search\"");
  coco_problem_t *problem;

  coco_reset_seeds();
  while ((problem = coco_suite_get_next_problem(suite, observer)) != NULL) {
    const size_t n = coco_problem_get_dimension(problem);
    const size_t m = coco_problem_get_number_of_objectives(problem);
    const size_t k = coco_problem_get_number_of_constraints(problem);
    const size_t budget = 50 * n;
    const double *lower = coco_problem_get_smallest_values_of_interest(problem);
    const double *upper = coco_problem_get_largest_values_of_interest(problem);
    double *x = (double *)malloc(n * sizeof(double));
    double *y = (double *)malloc(m * sizeof(double));
    double *c = (double *)malloc((k > 0 ? k : 1) * sizeof(double));
    size_t i, j;

    rng_state = coco_problem_get_suite_dep_index(problem);
    coco_problem_get_initial_solution(problem, x);
    coco_evaluate_function(problem, x, y);
    for (j = 0; j < budget; ++j) {
      if (j == budget / 2)
        coco_observer_signal_restart(observer, problem);
      for (i = 0; i < n; ++i)
        x[i] = lower[i] + rng_uniform() * (upper[i] - lower[i]);
      if (k > 0)
        coco_evaluate_constraint(problem, x, c);
      coco_evaluate_function(problem, x, y);
    }
    /* The loggers have to cope with non-finite values */
    x[0] = NAN;
    if (k > 0)
      coco_evaluate_constraint(problem, x, c);
    coco_evaluate_function(problem, x, y);
    x[0] = INFINITY;
    if (k > 0)
      coco_evaluate_constraint(problem, x, c);
    coco_evaluate_function(problem, x, y);
    coco_recommend_solution(problem, x);

    free(x);
    free(y);
    free(c);
  }
  coco_observer_free(observer);
  coco_suite_free(suite);
  return EXIT_SUCCESS;
}

int main(int argc, char **argv) {
  coco_set_log_level("warning");
  if (argc >= 3 && strcmp(argv[1], "values") == 0) {
    dump = argc >= 4 && strcmp(argv[3], "--dump") == 0;
    return values(argv[2]);
  } else if (argc == 5 && strcmp(argv[1], "observe") == 0) {
    return observe(argv[2], argv[3], argv[4]);
  }
  fprintf(stderr, "Usage: %s values <suite> [--dump]\n", argv[0]);
  fprintf(stderr, "       %s observe <suite> <observer> <suite options>\n", argv[0]);
  return EXIT_FAILURE;
}
