#include "coco.h"
#include "minunit.h"

/**
 * Tests the function coco_suite_get_problem.
 */
MU_TEST(test_coco_suite_get_problem) {

  coco_suite_t *suite = coco_suite("bbob-biobj", "instances: 1-10", "dimensions: 5");
  coco_problem_t *problem;

  problem = coco_suite_get_problem(suite, 0);
  mu_check(problem == NULL);

  problem = coco_suite_get_problem(suite, 1200);
  mu_check(problem != NULL);

  coco_problem_free(problem);
  coco_suite_free(suite);
}

/**
 * Tests the function coco_suite_encode_problem_index.
 */
MU_TEST(test_coco_suite_encode_problem_index) {

  coco_suite_t *suite;
  size_t index;

  suite = coco_suite("bbob", "year: 0000", NULL);
  index = coco_suite_encode_problem_index(suite, 13, 0, 10);
  mu_check(index == 205);
  coco_suite_free(suite);
}

/**
 * Tests the functions that map function, dimension and instance indices to their values.
 */
MU_TEST(test_coco_suite_get_from_index) {

  coco_suite_t *suite = coco_suite("bbob", "instances: 3,5,7", "");

  mu_check(coco_suite_get_function_from_function_index(suite, 0) == 1);
  mu_check(coco_suite_get_function_from_function_index(suite, 23) == 24);
  mu_check(coco_suite_get_dimension_from_dimension_index(suite, 0) == 2);
  mu_check(coco_suite_get_dimension_from_dimension_index(suite, 4) == 20);
  mu_check(coco_suite_get_instance_from_instance_index(suite, 0) == 3);
  mu_check(coco_suite_get_instance_from_instance_index(suite, 1) == 5);
  mu_check(coco_suite_get_instance_from_instance_index(suite, 2) == 7);

  coco_suite_free(suite);
}

/**
 * Run all tests in this file.
 */
int main(void) {
  MU_RUN_TEST(test_coco_suite_encode_problem_index);
  MU_RUN_TEST(test_coco_suite_get_problem);
  MU_RUN_TEST(test_coco_suite_get_from_index);

  MU_REPORT();

  return MU_EXIT_CODE;
}

