/**
 * Calls of the public C interface (coco.h) that must end in coco_error(), which exits the process.
 *
 * Each case is a separate ctest test that passes if the expected "COCO FATAL ERROR" message is printed (see
 * CMakeLists.txt). If the call returns, the program says so and the message is missing. The tests run in the
 * build directory, where the observers leave their output.
 */
#include <stdio.h>
#include <string.h>

#include "coco.h"

static int run(const char *name) {
  coco_set_log_level("error");

  if (strcmp(name, "unknown_suite") == 0) {
    coco_suite("no_such_suite", "", "");
  } else if (strcmp(name, "allocate_zero") == 0) {
    coco_allocate_memory(0);
  } else if (strcmp(name, "function_index") == 0) {
    coco_suite_t *suite = coco_suite("bbob", "", "");
    coco_suite_get_function_from_function_index(suite, 24);
  } else if (strcmp(name, "dimension_index") == 0) {
    coco_suite_t *suite = coco_suite("bbob", "", "");
    coco_suite_get_dimension_from_dimension_index(suite, 6);
  } else if (strcmp(name, "instance_index") == 0) {
    coco_suite_t *suite = coco_suite("bbob", "instances: 1-3", "");
    coco_suite_get_instance_from_instance_index(suite, 3);
  } else if (strcmp(name, "fvalues_of_interest") == 0) {
    coco_suite_t *suite = coco_suite("bbob", "", "");
    coco_problem_get_largest_fvalues_of_interest(coco_suite_get_problem(suite, 0));
  } else if (strcmp(name, "remove_other_observer") == 0) {
    coco_suite_t *suite = coco_suite("bbob", "instances: 1", "dimensions: 2");
    coco_observer_t *bbob = coco_observer("bbob", "outer_folder: test_coco_api_errors_exdata");
    coco_observer_t *toy = coco_observer("toy", "outer_folder: test_coco_api_errors_exdata");
    coco_problem_t *problem = coco_suite_get_next_problem(suite, bbob);
    coco_problem_remove_observer(problem, toy);
  } else if (strcmp(name, "feed_solution_order") == 0) {
    coco_suite_t *suite = coco_suite("bbob-biobj", "instances: 1", "dimensions: 2");
    coco_observer_t *observer =
        coco_observer("bbob-biobj", "outer_folder: test_coco_api_errors_exdata log_nondominated: read");
    coco_problem_t *problem = coco_suite_get_next_problem(suite, observer);
    const double y[2] = {1.0, 1.0};
    coco_logger_biobj_feed_solution(problem, 2, y);
    coco_logger_biobj_feed_solution(problem, 1, y);
  } else if (strcmp(name, "feed_solution_unobserved") == 0) {
    coco_suite_t *suite = coco_suite("bbob-biobj", "instances: 1", "dimensions: 2");
    coco_problem_t *problem = coco_suite_get_next_problem(suite, NULL);
    const double y[2] = {1.0, 1.0};
    coco_logger_biobj_feed_solution(problem, 1, y);
  } else {
    printf("unknown case %s\n", name);
    return 2;
  }
  printf("no error in case %s\n", name);
  return 1;
}

int main(int argc, char **argv) {
  if (argc != 2) {
    printf("usage: %s <case>\n", argv[0]);
    return 2;
  }
  return run(argv[1]);
}
