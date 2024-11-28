/**
 * @file observer_sbox_cost.cpp
 * @brief Implementation of the sbox_cost observer.
 */

#include "coco.h"
#include "coco_utilities.cpp"

// static coco_problem_t* logger_sbox_cost(coco_observer_t *observer, coco_problem_t* problem);
// static void logger_sbox_cost_free(void* logger);

/**
 * @brief The sbox_cost observer data type.
 */
typedef struct {
  /* TODO: Can be used to store variables that need to be accessible during one run (i.e. for multiple
   * problems). For example, the following global variables from logger_sbox_cost.c could be stored here: */
  size_t current_dim;
  size_t current_fun_id;
  /* ... and others */
} observer_sbox_cost_data_t;

/**
 * @brief Initializes the sbox_cost observer.
 */
static void observer_sbox_cost(coco_observer_t *observer, char const* options, coco_option_keys_t** option_keys) {
  // FIXME: Use sbox_cost logger once ported to new logging code
  observer->logger_allocate_function = logger_bbob;
  observer->logger_free_function = logger_bbob_free;
  observer->data_free_function = nullptr;
  observer->data = nullptr;

  *option_keys = nullptr;

  (void) options; /* To silence the compiler */
}
