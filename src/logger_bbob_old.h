#pragma once
#include "coco.h"
#include "coco_internal.h"
#include "coco_observer.h"

#include <stdio.h>

/**
 * @brief The bbob_old logger data type.
 */
typedef struct {
  coco_observer_t *observer;
  int is_initialized;
  /*char *path;// relative path to the data folder. //Wassim: now fetched from the observer */
  /*const char *alg_name; the alg name, for now, temporarily the same as the path. Wassim: Now in the observer */
  FILE *index_file; /* index file */
  FILE *fdata_file; /* function value aligned data file */
  FILE *tdata_file; /* number of function evaluations aligned data file */
  FILE *rdata_file; /* restart info data file */
  size_t number_of_evaluations;
  size_t number_of_evaluations_constraints;
  double best_fvalue;
  double last_fvalue;
  short written_last_eval; /* allows writing the data of the final fun eval in the .tdat file if not already written by
                              the t_trigger*/
  double *best_solution;
  /* The following are to only pass data as a parameter in the free function. The
   * interface should probably be the same for all free functions so passing the
   * problem as a second parameter is not an option even though we need info
   * form it.*/
  size_t function_id; /*TODO: consider changing name*/
  size_t instance_id;
  size_t number_of_variables;
  size_t number_of_integer_variables;
  int log_discrete_as_int; /**< @brief Whether to output discrete variables in int or double format. */
  double optimal_fvalue;
  char *suite_name;

  coco_observer_log_targets_t *targets;     /**< @brief Triggers based on target values. */
  coco_observer_evaluations_t *evaluations; /**< @brief Triggers based on the number of evaluations. */
} logger_bbob_old_data_t;

coco_problem_t *logger_bbob_old(coco_observer_t *observer, coco_problem_t *inner_problem);
void logger_bbob_old_free(void *stuff);
void logger_bbob_old_evaluate(coco_problem_t *problem, const double *x, double *y);
void logger_bbob_old_initialize(logger_bbob_old_data_t *logger, coco_problem_t *inner_problem);
void logger_bbob_old_openIndexFile(logger_bbob_old_data_t *logger, const char *folder_path,
                                   const char *indexFile_prefix, const char *function_id, const char *dataFile_path,
                                   const char *suite_name);
