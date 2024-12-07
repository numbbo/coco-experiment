#pragma once
#include "coco.h"
#include "coco_internal.h"

/**
 * @brief The bbob observer data type.
 *
 * There is a cyclic reference between the bbob_logger and the bbob_observer (through the observer's
 * observed_problem and the logger's data, which points to the observer). This is needed to be able
 * to free both objects without problems.
 */
typedef struct {
  coco_problem_t *observed_problem; /**< @brief Pointer to the observed problem (NULL if none is observed) */
  char *prefix;                     /**< @brief Prefix in the name of the info and data files */

  /* Store the information on the last function, dimension and dat file used, to be able to tell when the
   * logger needs to write in a new dat file.  */
  size_t last_function;  /**< @brief The function that was logged last */
  size_t last_dimension; /**< @brief The dimension that was logged last */
  char *last_dat_file;   /**< @brief The name of the .dat file that was last used for logging */
} observer_bbob_data_t;

void observer_bbob(coco_observer_t *observer, const char *options, coco_option_keys_t **option_keys);
void observer_bbob_data_free(void *stuff);
void logger_bbob_data_nullify_observer(void *logger_data);
void logger_bbob_signal_restart(coco_problem_t *problem);
void logger_bbob_free(void *logger);
coco_problem_t *logger_bbob(coco_observer_t *observer, coco_problem_t *problem);
