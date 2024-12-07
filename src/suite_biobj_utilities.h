#pragma once
#include "coco.h"
#include "coco_internal.h"

/**
 * @brief A structure containing information about the new instances.
 */
typedef struct {

  size_t **new_instances; /**< @brief A matrix of new instances (equal in form to suite_biobj_instances)
                                that needs to be used only when an instance that is not in
                                suite_biobj_instances is being invoked. */

  size_t max_new_instances; /**< @brief The maximal number of new instances. */

} suite_biobj_new_inst_t;

void suite_biobj_get_best_hyp_value(const int known_optima, const char *key, double *value);
coco_problem_t *coco_get_biobj_problem(const size_t function, const size_t dimension, const size_t instance,
                                       const coco_get_problem_function_t coco_get_problem_function,
                                       suite_biobj_new_inst_t **new_inst_data, const size_t num_new_instances,
                                       const size_t *dimensions, const size_t num_dimensions);
size_t suite_biobj_get_new_instance(suite_biobj_new_inst_t *new_inst_data, const size_t instance,
                                    const size_t instance1, const size_t *bbob_functions,
                                    const size_t num_bbob_functions, const size_t *sel_bbob_functions,
                                    const size_t num_sel_bbob_functions, const size_t *dimensions,
                                    const size_t num_dimensions);
int suite_biobj_check_inst_consistency(const size_t dimension, size_t function1, size_t instance1, size_t function2,
                                       size_t instance2);
void suite_biobj_new_inst_free(void *stuff);
extern const size_t suite_biobj_instances[][3];
