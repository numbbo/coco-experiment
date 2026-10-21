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
 * Tests that coco_suite_decode_problem_index sets all indices to 0 for an index that is too large.
 */
MU_TEST(test_coco_suite_decode_problem_index) {

  coco_suite_t *suite = coco_suite("bbob", "year: 0000", NULL);
  size_t function_idx = 7, dimension_idx = 7, instance_idx = 7;

  coco_suite_decode_problem_index(suite, 205, &function_idx, &dimension_idx, &instance_idx);
  mu_check((function_idx == 13) && (dimension_idx == 0) && (instance_idx == 10));

  coco_set_log_level("error");
  coco_suite_decode_problem_index(suite, 1000000, &function_idx, &dimension_idx, &instance_idx);
  coco_set_log_level("info");
  mu_check((function_idx == 0) && (dimension_idx == 0) && (instance_idx == 0));

  coco_suite_free(suite);
}

/**
 * Tests that a suite left without dimensions by the filters is not returned.
 */
MU_TEST(test_coco_suite_empty_after_filtering) {

  coco_set_log_level("error");
  mu_check(coco_suite("bbob", "", "dimensions: 7") == NULL);
  coco_set_log_level("info");
}

/**
 * Run all tests in this file.
 */
int main(void) {
  MU_RUN_TEST(test_coco_suite_encode_problem_index);
  MU_RUN_TEST(test_coco_suite_get_problem);
  MU_RUN_TEST(test_coco_suite_get_from_index);
  MU_RUN_TEST(test_coco_suite_decode_problem_index);
  MU_RUN_TEST(test_coco_suite_empty_after_filtering);

  MU_REPORT();

  return MU_EXIT_CODE;
}

