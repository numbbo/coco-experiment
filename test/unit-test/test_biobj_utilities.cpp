#include "minunit.h"

#include "coco.cpp"

/**
 * Tests the creation of new instances of the bi-objective suites.
 */
MU_TEST(test_biobj_utilities_instances) {

  coco_suite_t *suite = coco_suite("bbob-biobj", "instances: 16-20", "");
  coco_problem_t *problem;

  problem = coco_suite_get_next_problem(suite, NULL);
  mu_check(ProblemGetSuiteDepInstance(FromC(problem)) == 16);

  problem = coco_suite_get_next_problem(suite, NULL);
  mu_check(ProblemGetSuiteDepInstance(FromC(problem)) == 17);

  problem = coco_suite_get_next_problem(suite, NULL);
  mu_check(ProblemGetSuiteDepInstance(FromC(problem)) == 18);

  problem = coco_suite_get_next_problem(suite, NULL);
  mu_check(ProblemGetSuiteDepInstance(FromC(problem)) == 19);

  problem = coco_suite_get_next_problem(suite, NULL);
  mu_check(ProblemGetSuiteDepInstance(FromC(problem)) == 20);

  coco_suite_free(suite);
}

/**
 * Run all tests in this file.
 */
int main(void) {
  MU_RUN_TEST(test_biobj_utilities_instances);

  MU_REPORT();

  return MU_EXIT_CODE;
}
