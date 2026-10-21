#include "minunit.h"

#include "coco.cpp"

/**
 * Tests the creation of new instances of the bi-objective suites.
 */
MU_TEST(test_biobj_utilities_instances) {

  coco_suite_t *suite = coco_suite("bbob-biobj", "instances: 16-20", "");
  coco_problem_t *problem;

  problem = coco_suite_get_next_problem(suite, NULL);
  mu_check(FromC(problem)->suite_dep_instance() == 16);

  problem = coco_suite_get_next_problem(suite, NULL);
  mu_check(FromC(problem)->suite_dep_instance() == 17);

  problem = coco_suite_get_next_problem(suite, NULL);
  mu_check(FromC(problem)->suite_dep_instance() == 18);

  problem = coco_suite_get_next_problem(suite, NULL);
  mu_check(FromC(problem)->suite_dep_instance() == 19);

  problem = coco_suite_get_next_problem(suite, NULL);
  mu_check(FromC(problem)->suite_dep_instance() == 20);

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
