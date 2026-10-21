#include "minunit.h"

#include "coco.cpp"

/** The tolerance for the objective value at the best parameter, which may be computed on a different path. */
static const double precision = 1e-8;

/**
 * Checks that the best parameter of every problem of the suite is an optimum: it lies in the region of interest
 * and its objective value (without noise) is the best value. Prints the problems that fail and returns their
 * number.
 */
static size_t check_best_parameters(const char *suite_name, const char *suite_options) {
  coco_suite_t *suite = coco_suite(suite_name, "", suite_options);
  coco_problem_t *problem;
  size_t number_of_problems = 0, number_of_failures = 0;

  while ((problem = coco_suite_get_next_problem(suite, nullptr)) != nullptr) {
    const size_t dimension = coco_problem_get_dimension(problem);
    const double *lower = coco_problem_get_smallest_values_of_interest(problem);
    const double *upper = coco_problem_get_largest_values_of_interest(problem);
    std::vector<double> x(dimension);
    double y;
    bool in_region_of_interest = true;

    ++number_of_problems;
    if (!coco_problem_get_best_parameter(problem, x.data())) {
      printf("%s: no best parameter\n", coco_problem_get_id(problem));
      ++number_of_failures;
      continue;
    }
    for (size_t i = 0; i < dimension; ++i)
      in_region_of_interest = in_region_of_interest && lower[i] <= x[i] && x[i] <= upper[i];

    coco_evaluate_function(problem, x.data(), &y);
    if (problem->is_noisy())
      y = problem->get_last_noise_free_values()[0];
    const double best_value = coco_problem_get_best_value(problem);

    if (!in_region_of_interest || !(fabs(y - best_value) <= precision)) {
      printf("%s: f(best parameter) = %.17g, best value = %.17g%s\n", coco_problem_get_id(problem), y, best_value,
             in_region_of_interest ? "" : ", best parameter outside the region of interest");
      ++number_of_failures;
    }
  }
  coco_suite_free(suite);
  /* An empty suite would pass without testing anything */
  return number_of_problems > 0 ? number_of_failures : 1;
}

MU_TEST(test_bbob) { mu_check(check_best_parameters("bbob", "") == 0); }

MU_TEST(test_bbob_boxed) { mu_check(check_best_parameters("bbob-boxed", "") == 0); }

MU_TEST(test_bbob_noisy) { mu_check(check_best_parameters("bbob-noisy", "") == 0); }

MU_TEST(test_bbob_largescale) { mu_check(check_best_parameters("bbob-largescale", "dimensions: 20,40,80") == 0); }

MU_TEST(test_bbob_mixint) { mu_check(check_best_parameters("bbob-mixint", "") == 0); }

MU_TEST(test_bbob_constrained) {
  mu_check(check_best_parameters("bbob-constrained", "") == 0);
  mu_check(check_best_parameters("bbob-constrained-active-only", "dimensions: 2,10") == 0);
  mu_check(check_best_parameters("bbob-constrained-no-disguise", "dimensions: 2,10") == 0);
}

MU_TEST(test_toy) { mu_check(check_best_parameters("toy", "") == 0); }

/**
 * Multi-objective problems have no best parameter.
 */
MU_TEST(test_biobj) {
  const char *suite_names[] = {"bbob-biobj", "bbob-biobj-ext", "bbob-biobj-mixint"};
  for (const char *suite_name : suite_names) {
    coco_suite_t *suite = coco_suite(suite_name, "", "dimensions: 2 instance_indices: 1");
    coco_problem_t *problem;
    size_t number_of_problems = 0, number_with_best_parameter = 0;
    while ((problem = coco_suite_get_next_problem(suite, nullptr)) != nullptr) {
      std::vector<double> x(coco_problem_get_dimension(problem));
      ++number_of_problems;
      number_with_best_parameter += coco_problem_get_best_parameter(problem, x.data());
    }
    coco_suite_free(suite);
    mu_check(number_of_problems > 0);
    mu_check(number_with_best_parameter == 0);
  }
}

/**
 * Run all tests in this file.
 */
int main(void) {
  coco_set_log_level("error");

  MU_RUN_TEST(test_bbob);
  MU_RUN_TEST(test_bbob_boxed);
  MU_RUN_TEST(test_bbob_noisy);
  MU_RUN_TEST(test_bbob_largescale);
  MU_RUN_TEST(test_bbob_mixint);
  MU_RUN_TEST(test_bbob_constrained);
  MU_RUN_TEST(test_toy);
  MU_RUN_TEST(test_biobj);

  MU_REPORT();

  return MU_EXIT_CODE;
}
