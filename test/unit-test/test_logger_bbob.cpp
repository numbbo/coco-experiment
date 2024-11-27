#include <stdio.h>

#include "minunit.h"
#include "about_equal.h"

#include "coco.cpp"

/**
 * Counts the number of lines in a file.
 */
int count_lines_in_file(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        coco_error("Error opening file %s", filename);
        return -1; 
    }
    int lines = 0;
    int ch;
    while ((ch = fgetc(file)) != EOF) {
        if (ch == '\n') {
            lines++;
        }
    }
    fclose(file);
    return lines;
}

/**
 * Tests the logarithmic and linear triggers via the bbob logger.
 */
MU_TEST(test_logger_bbob_triggers) {

  coco_suite_t *suite;
  coco_observer_t *observer;
  coco_problem_t *problem;

  double *x = coco_allocate_vector(2);
  double *y = coco_allocate_vector(1);
  double target;

  logger_bbob_data_t *logger;

  /* Using only the logarithmic performance targets (the linear target precision needs to be high) */
  suite = coco_suite("bbob", "", "dimensions: 2 function_indices: 1 instance_indices: 1");
  observer = coco_observer("bbob", "number_target_triggers: 1 lin_target_precision: 100");
  /* Use the 2-D sphere function */
  problem = coco_suite_get_next_problem(suite, observer);
  logger = (logger_bbob_data_t *) coco_problem_transformed_get_data(problem);

  x[0] = 0; x[1] = 0;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 80.88209408));
  coco_observer_targets_trigger(logger->targets, logger->best_found_value - logger->optimal_value);
  mu_check(about_equal_value(logger->best_found_value - logger->optimal_value, 1.40209408));
  target = coco_observer_targets_get_last_target((coco_observer_targets_t *)logger->targets);
  mu_check(about_equal_value(target, 10));

  x[0] = 0; x[1] = -0.5;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.97529408));
  coco_observer_targets_trigger(logger->targets, logger->best_found_value - logger->optimal_value);
  mu_check(about_equal_value(logger->best_found_value - logger->optimal_value, 0.49529408));
  target = coco_observer_targets_get_last_target(logger->targets);
  mu_check(about_equal_value(target, 1));

  x[0] = 0; x[1] = -1;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.56849408));
  coco_observer_targets_trigger(logger->targets, logger->best_found_value - logger->optimal_value);
  mu_check(about_equal_value(logger->best_found_value - logger->optimal_value, 0.08849408));
  target = coco_observer_targets_get_last_target(logger->targets);
  mu_check(about_equal_value(target, 0.1));

  x[0] = 0; x[1] = -1.1;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.54713408));
  coco_observer_targets_trigger(logger->targets, logger->best_found_value - logger->optimal_value);
  mu_check(about_equal_value(logger->best_found_value - logger->optimal_value, 0.06713408));
  target = coco_observer_targets_get_last_target(logger->targets);
  mu_check(about_equal_value(target, 0.1));

  x[0] = 0.2; x[1] = -1.1;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.48601408));
  coco_observer_targets_trigger(logger->targets, logger->best_found_value - logger->optimal_value);
  mu_check(about_equal_value(logger->best_found_value - logger->optimal_value, 0.00601408));
  target = coco_observer_targets_get_last_target(logger->targets);
  mu_check(about_equal_value(target, 0.01));

  coco_observer_free(observer);
  coco_suite_free(suite);

  /* Using the linear performance targets */
  suite = coco_suite("bbob", "", "dimensions: 2 function_indices: 1 instance_indices: 1");
  suite->known_optima = 0;
  observer = coco_observer("bbob", "number_target_triggers: 1 lin_target_precision: 1e-3 log_target_precision: 1e-5");
  /* Use the 2-D sphere function */
  problem = coco_suite_get_next_problem(suite, observer);
  logger = (logger_bbob_data_t *) coco_problem_transformed_get_data(problem);

  x[0] = 0; x[1] = 0;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 80.88209408));
  coco_observer_targets_trigger(logger->targets, logger->best_found_value - logger->optimal_value);
  mu_check(about_equal_value(logger->best_found_value - logger->optimal_value, 1.40209408));
  target = coco_observer_targets_get_last_target((coco_observer_targets_t *)logger->targets);
  mu_check(about_equal_value(target, 1.403));

  x[0] = 0; x[1] = -0.5;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.97529408));
  coco_observer_targets_trigger(logger->targets, logger->best_found_value - logger->optimal_value);
  mu_check(about_equal_value(logger->best_found_value - logger->optimal_value, 0.49529408));
  target = coco_observer_targets_get_last_target(logger->targets);
  mu_check(about_equal_value(target, 0.496));

  x[0] = 0; x[1] = -1;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.56849408));
  coco_observer_targets_trigger(logger->targets, logger->best_found_value - logger->optimal_value);
  mu_check(about_equal_value(logger->best_found_value - logger->optimal_value, 0.08849408));
  target = coco_observer_targets_get_last_target(logger->targets);
  mu_check(about_equal_value(target, 0.089));

  x[0] = 0; x[1] = -1.1;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.54713408));
  coco_observer_targets_trigger(logger->targets, logger->best_found_value - logger->optimal_value);
  mu_check(about_equal_value(logger->best_found_value - logger->optimal_value, 0.06713408));
  target = coco_observer_targets_get_last_target(logger->targets);
  mu_check(about_equal_value(target, 0.068));

  x[0] = 0.2; x[1] = -1.1;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.48601408));
  coco_observer_targets_trigger(logger->targets, logger->best_found_value - logger->optimal_value);
  mu_check(about_equal_value(logger->best_found_value - logger->optimal_value, 0.00601408));
  target = coco_observer_targets_get_last_target(logger->targets);
  mu_check(about_equal_value(target, 0.007));

  coco_observer_free(observer);
  coco_suite_free(suite);

  /* Testing the targets that fall below log_target_precision */
  suite = coco_suite("bbob", "", "dimensions: 2 function_indices: 1 instance_indices: 1");
  suite->known_optima = 1;
  observer = coco_observer("bbob", "number_target_triggers: 1 lin_target_precision: 10 log_target_precision: 1e-5");
  /* Use the 2-D sphere function */
  problem = coco_suite_get_next_problem(suite, observer);
  logger = (logger_bbob_data_t *) coco_problem_transformed_get_data(problem);

  x[0] = 0.252; x[1] = -1.156;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.48000128));
  coco_observer_targets_trigger(logger->targets, logger->best_found_value - logger->optimal_value);
  mu_check(about_equal_value(logger->best_found_value - logger->optimal_value, 1.28 * 1e-6));
  target = coco_observer_targets_get_last_target(logger->targets);
  mu_check(about_equal_value(target, 1e-5));

  coco_observer_free(observer);
  coco_suite_free(suite);

  coco_free_memory(x);
  coco_free_memory(y);
}

/**
 * Test whether the last line is correctly printed to the various log files. 
 */
MU_TEST(test_logger_bbob_last_line) {

  coco_suite_t *suite;
  coco_observer_t *observer;
  coco_problem_t *problem;

  double *x = coco_allocate_vector(2);
  double *y = coco_allocate_vector(1);

  int i = 0;
  char *filename, *folder;

  suite = coco_suite("bbob", "", "dimensions: 2 function_indices: 1 instance_indices: 1");
  observer = coco_observer("bbob", "number_evaluation_triggers: 4 base_evaluation_triggers: 100");
  /* Use the 2-D sphere function */
  problem = coco_suite_get_next_problem(suite, observer);
  folder = coco_strdup(observer->result_folder);

  x[0] = 0; x[1] = 0;
  coco_evaluate_function(problem, x, y);
  coco_recommend_solution(problem, x);
  mu_check(about_equal_value(y[0], 80.88209408));

  x[0] = 0; x[1] = -0.5;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.97529408));

  x[0] = 0.2; x[1] = -1.1;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.48601408));

  for (i = 0; i < 5; i++) {
    x[0] = 0; x[1] = 0;
    coco_evaluate_function(problem, x, y);
  }

  coco_observer_free(observer);
  coco_suite_free(suite);

  filename =
      coco_strdupf("%s%s%s%s%s", folder, coco_path_separator, "data_f1", coco_path_separator, "bbobexp_f1_DIM2.dat");
  mu_check(count_lines_in_file(filename) == 5);
  coco_free_memory(filename);

  filename = coco_strdupf("%s%s%s%s%s", folder, coco_path_separator, "data_f1", 
    coco_path_separator, "bbobexp_f1_DIM2.tdat");
  mu_check(count_lines_in_file(filename) == 9);
  coco_free_memory(filename);

  filename = coco_strdupf("%s%s%s%s%s", folder, coco_path_separator, "data_f1", 
    coco_path_separator, "bbobexp_f1_DIM2.mdat");
  mu_check(count_lines_in_file(filename) == 3);

  coco_free_memory(x);
  coco_free_memory(y);
  coco_free_memory(folder);
  coco_free_memory(filename);
}


/**
 * Run all tests in this file.
 */
int main(void) {
  MU_RUN_TEST(test_logger_bbob_triggers);
  MU_RUN_TEST(test_logger_bbob_last_line);

  MU_REPORT();

  return MU_EXIT_CODE;
}
