#pragma once
#include "coco.h"
#include "observer_bbob.h"

typedef struct logger_bbob_data_s logger_bbob_data_t;

coco_problem_t *logger_bbob(coco_observer_t *observer, coco_problem_t *inner_problem);
void logger_bbob_signal_restart(coco_problem_t *problem);
void logger_bbob_data_nullify_observer(void *stuff);
void logger_bbob_free(void *stuff);
void logger_bbob_recommend(coco_problem_t *problem, const double *x);
void logger_bbob_evaluate(coco_problem_t *problem, const double *x, double *y);
void logger_bbob_initialize(logger_bbob_data_t *logger, int is_opt_known);
void logger_bbob_open_info_file(logger_bbob_data_t *logger, const char *folder, const char *function_string,
                                const char *data_file_name, const char *suite_name, int start_new_line);
/*
void logger_bbob_open_data_file(FILE **data_file, const char *path, const char *file_name, const char *file_extension);
void logger_bbob_open_file(FILE **file, const char *file_path);
void logger_bbob_output(FILE *data_file, logger_bbob_data_t *logger, const double *x, double current_value,
                        double best_value, const double *constraints);
                        */
int logger_bbob_start_new_line(coco_observer_t *observer, size_t current_dimension, size_t current_function);
int logger_bbob_single_digit_constraint(const double c);
extern const char *logger_bbob_header;
extern const char *logger_bbob_data_format;
