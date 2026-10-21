#include <stdio.h>

#include "minunit.h"
#include "about_equal.h"

#include "coco.cpp"

namespace coco {

/**
 * Counts the number of lines in a file.
 */
int count_lines_in_file(const std::string &filename) {
    FILE *file = fopen(filename.c_str(), "r");
    if (file == NULL) {
        coco_error("Error opening file %s", filename.c_str());
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

  LoggerBbob *logger;

  /* Using only the logarithmic performance targets (the linear target precision needs to be high) */
  suite = coco_suite("bbob", "", "dimensions: 2 function_indices: 1 instance_indices: 1");
  observer = coco_observer("bbob", "number_target_triggers: 1 lin_target_precision: 100");
  /* Use the 2-D sphere function */
  problem = coco_suite_get_next_problem(suite, observer);
  logger = dynamic_cast<LoggerBbob *>(FromC(problem));
  mu_check(logger != nullptr);

  x[0] = 0; x[1] = 0;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 80.88209408));
  logger->targets().Trigger(logger->best_found_value() - logger->optimal_value());
  mu_check(about_equal_value(logger->best_found_value() - logger->optimal_value(), 1.40209408));
  target = logger->targets().LastTarget();
  mu_check(about_equal_value(target, 10));

  x[0] = 0; x[1] = -0.5;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.97529408));
  logger->targets().Trigger(logger->best_found_value() - logger->optimal_value());
  mu_check(about_equal_value(logger->best_found_value() - logger->optimal_value(), 0.49529408));
  target = logger->targets().LastTarget();
  mu_check(about_equal_value(target, 1));

  x[0] = 0; x[1] = -1;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.56849408));
  logger->targets().Trigger(logger->best_found_value() - logger->optimal_value());
  mu_check(about_equal_value(logger->best_found_value() - logger->optimal_value(), 0.08849408));
  target = logger->targets().LastTarget();
  mu_check(about_equal_value(target, 0.1));

  x[0] = 0; x[1] = -1.1;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.54713408));
  logger->targets().Trigger(logger->best_found_value() - logger->optimal_value());
  mu_check(about_equal_value(logger->best_found_value() - logger->optimal_value(), 0.06713408));
  target = logger->targets().LastTarget();
  mu_check(about_equal_value(target, 0.1));

  x[0] = 0.2; x[1] = -1.1;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.48601408));
  logger->targets().Trigger(logger->best_found_value() - logger->optimal_value());
  mu_check(about_equal_value(logger->best_found_value() - logger->optimal_value(), 0.00601408));
  target = logger->targets().LastTarget();
  mu_check(about_equal_value(target, 0.01));

  coco_observer_free(observer);
  coco_suite_free(suite);

  /* Using the linear performance targets */
  suite = coco_suite("bbob", "", "dimensions: 2 function_indices: 1 instance_indices: 1");
  FromC(suite)->set_known_optima(false);
  observer = coco_observer("bbob", "number_target_triggers: 1 lin_target_precision: 1e-3 log_target_precision: 1e-5");
  /* Use the 2-D sphere function */
  problem = coco_suite_get_next_problem(suite, observer);
  logger = dynamic_cast<LoggerBbob *>(FromC(problem));
  mu_check(logger != nullptr);

  x[0] = 0; x[1] = 0;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 80.88209408));
  logger->targets().Trigger(logger->best_found_value() - logger->optimal_value());
  mu_check(about_equal_value(logger->best_found_value() - logger->optimal_value(), 1.40209408));
  target = logger->targets().LastTarget();
  mu_check(about_equal_value(target, 1.403));

  x[0] = 0; x[1] = -0.5;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.97529408));
  logger->targets().Trigger(logger->best_found_value() - logger->optimal_value());
  mu_check(about_equal_value(logger->best_found_value() - logger->optimal_value(), 0.49529408));
  target = logger->targets().LastTarget();
  mu_check(about_equal_value(target, 0.496));

  x[0] = 0; x[1] = -1;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.56849408));
  logger->targets().Trigger(logger->best_found_value() - logger->optimal_value());
  mu_check(about_equal_value(logger->best_found_value() - logger->optimal_value(), 0.08849408));
  target = logger->targets().LastTarget();
  mu_check(about_equal_value(target, 0.089));

  x[0] = 0; x[1] = -1.1;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.54713408));
  logger->targets().Trigger(logger->best_found_value() - logger->optimal_value());
  mu_check(about_equal_value(logger->best_found_value() - logger->optimal_value(), 0.06713408));
  target = logger->targets().LastTarget();
  mu_check(about_equal_value(target, 0.068));

  x[0] = 0.2; x[1] = -1.1;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.48601408));
  logger->targets().Trigger(logger->best_found_value() - logger->optimal_value());
  mu_check(about_equal_value(logger->best_found_value() - logger->optimal_value(), 0.00601408));
  target = logger->targets().LastTarget();
  mu_check(about_equal_value(target, 0.007));

  coco_observer_free(observer);
  coco_suite_free(suite);

  /* Testing the targets that fall below log_target_precision */
  suite = coco_suite("bbob", "", "dimensions: 2 function_indices: 1 instance_indices: 1");
  FromC(suite)->set_known_optima(true);
  observer = coco_observer("bbob", "number_target_triggers: 1 lin_target_precision: 10 log_target_precision: 1e-5");
  /* Use the 2-D sphere function */
  problem = coco_suite_get_next_problem(suite, observer);
  logger = dynamic_cast<LoggerBbob *>(FromC(problem));
  mu_check(logger != nullptr);

  x[0] = 0.252; x[1] = -1.156;
  coco_evaluate_function(problem, x, y);
  mu_check(about_equal_value(y[0], 79.48000128));
  logger->targets().Trigger(logger->best_found_value() - logger->optimal_value());
  mu_check(about_equal_value(logger->best_found_value() - logger->optimal_value(), 1.28 * 1e-6));
  target = logger->targets().LastTarget();
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
  std::string folder;

  coco_set_log_level("debug");

  suite = coco_suite("bbob", "", "dimensions: 2 function_indices: 1 instance_indices: 1");
  observer = coco_observer("bbob", "number_target_triggers: 100 number_evaluation_triggers: 20 base_evaluation_triggers: 100");
  /* Use the 2-D sphere function */
  problem = coco_suite_get_next_problem(suite, observer);
  folder = FromC(observer)->GetResultFolder();

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

  mu_check(count_lines_in_file(JoinPath(folder, "data_f1", "bbobexp_f1_DIM2.dat")) == 5);

  mu_check(count_lines_in_file(JoinPath(folder, "data_f1", "bbobexp_f1_DIM2.tdat")) == 9);

  mu_check(count_lines_in_file(JoinPath(folder, "data_f1", "bbobexp_f1_DIM2.mdat")) == 3);

  /* Repeat the check with different trigger values */
  suite = coco_suite("bbob", "", "dimensions: 2 function_indices: 1 instance_indices: 1");
  observer = coco_observer("bbob", "number_target_triggers: 20 number_evaluation_triggers: 4 base_evaluation_triggers: 100");
  /* Use the 2-D sphere function */
  problem = coco_suite_get_next_problem(suite, observer);
  folder = FromC(observer)->GetResultFolder();

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

  mu_check(count_lines_in_file(JoinPath(folder, "data_f1", "bbobexp_f1_DIM2.dat")) == 5);

  mu_check(count_lines_in_file(JoinPath(folder, "data_f1", "bbobexp_f1_DIM2.tdat")) == 5);

  mu_check(count_lines_in_file(JoinPath(folder, "data_f1", "bbobexp_f1_DIM2.mdat")) == 3);

  coco_free_memory(x);
  coco_free_memory(y);
}

/**
 * Run all tests in this file.
 */
/**
 * Tests that the bbob-old logger can start logging after its observer has been freed.
 */
MU_TEST(test_logger_bbob_old_evaluate_after_observer_free) {

  coco_suite_t *suite = coco_suite("bbob", "instances: 1", "dimensions: 2 function_indices: 1");
  coco_observer_t *observer = coco_observer("bbob-old", "result_folder: bbob_old_after_free");
  coco_problem_t *problem = coco_suite_get_next_problem(suite, observer);
  double x[2] = {0.0, 0.0};
  double y[1];

  mu_check(problem != NULL);
  coco_observer_free(observer);
  coco_evaluate_function(problem, x, y); /* The logger creates its files only now */
  coco_suite_free(suite);
}

} // namespace coco

int main(void) {
  MU_RUN_TEST(coco::test_logger_bbob_triggers);
  MU_RUN_TEST(coco::test_logger_bbob_last_line);
  MU_RUN_TEST(coco::test_logger_bbob_old_evaluate_after_observer_free);

  MU_REPORT();

  return MU_EXIT_CODE;
}
