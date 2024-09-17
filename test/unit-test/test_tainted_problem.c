/**
 * Tests different combinations of suite options
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "coco.h"

int main(int argc, const char **argv) {
  double x[100]; // Needlesly large
  double y[1];
  coco_problem_t *problem;

  if (argc != 2) {
    fprintf(stderr, "USAGE: %s <suite>", argv[0]);
    // Important! Must return 0 because non-zero return is considered a test success!
    return 0; 
  }
  coco_suite_t *suite = coco_suite(argv[1], NULL, "");
  coco_observer_t *observer = coco_observer("bbob", "");
  
  while ((problem = coco_suite_get_next_problem(suite, observer)) != NULL) {
    coco_problem_get_best_parameter(problem, x);
    coco_evaluate_function(problem, x, y);
  }

  coco_observer_free(observer);
  coco_suite_free(suite);

  return 0;
}
