#pragma once
#include "coco.h"

typedef struct logger_biobj_data_s logger_biobj_data_t;

coco_problem_t *logger_biobj(coco_observer_t *observer, coco_problem_t *inner_problem);
void logger_biobj_signal_restart(coco_problem_t *problem);
void logger_biobj_data_nullify_observer(void *stuff);
void logger_biobj_free(void *stuff);
void logger_biobj_finalize(logger_biobj_data_t *logger);
int coco_logger_biobj_feed_solution(coco_problem_t *problem, const size_t evaluation, const double *y);
void logger_biobj_recommend(coco_problem_t *problem, const double *x);
void logger_biobj_evaluate(coco_problem_t *problem, const double *x, double *y);
void logger_biobj_indicator_free(void *stuff);
/*
void logger_biobj_indicator_initialize_file(const logger_biobj_data_t *logger, const coco_observer_t *observer,
                                            const coco_problem_t *problem, const logger_biobj_indicator_t *indicator,
                                            FILE **f, const char *file_ending, const int output_targets);
int logger_biobj_tree_update(logger_biobj_data_t *logger, logger_biobj_avl_item_t *node_item);
size_t logger_biobj_tree_output(FILE *file, const avl_tree_t *tree, const size_t dim, const size_t num_int_vars,
                                const size_t num_obj, const int log_vars, const int precision_x, const int precision_f,
                                const int log_discrete_as_int);
int avl_tree_compare_by_eval_number(const logger_biobj_avl_item_t *item1, const logger_biobj_avl_item_t *item2,
                                    void *userdata);
int avl_tree_compare_by_last_objective(const logger_biobj_avl_item_t *item1, const logger_biobj_avl_item_t *item2,
                                       void *userdata);
void logger_biobj_node_free(logger_biobj_avl_item_t *item, void *userdata);
logger_biobj_avl_item_t *logger_biobj_node_create(const coco_problem_t *problem, const double *x, const double *y,
                                                  const double *constraints, const size_t evaluation_number,
                                                  const size_t dim, const size_t num_obj, const size_t num_const);
                                                  */
