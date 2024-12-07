#pragma once
#include "coco.h"

/**
 * @brief The bbob-biobj-ext suite data type.
 */
typedef struct {

  size_t **new_instances; /**< @brief A matrix of new instances (equal in form to suite_biobj_ext_instances)
                                that needs to be used only when an instance that is not in
                                suite_biobj_ext_instances is being invoked. */

  size_t max_new_instances; /**< @brief The maximal number of new instances. */

} suite_biobj_ext_t;

int check_consistency_of_instances(const size_t dimension, size_t function1, size_t instance1, size_t function2,
                                   size_t instance2);
coco_problem_t *suite_biobj_ext_get_problem(coco_suite_t *suite, const size_t function_idx, const size_t dimension_idx,
                                            const size_t instance_idx);
const char *suite_biobj_ext_get_instances_by_year(const int year);
coco_suite_t *suite_biobj_ext_initialize(void);
size_t suite_biobj_ext_get_new_instance(coco_suite_t *suite, const size_t instance, const size_t instance1,
                                        const size_t num_bbob_functions, const size_t *bbob_functions);
size_t suite_biobj_ext_get_new_instance(coco_suite_t *suite, const size_t instance, const size_t instance1,
                                        const size_t num_bbob_functions, const size_t *bbob_functions);
void suite_biobj_ext_free(void *stuff);
void suite_biobj_ext_free(void *stuff);
coco_suite_t *coco_suite_allocate(const char *suite_name, const size_t number_of_functions,
                                  const size_t number_of_dimensions, const size_t *dimensions,
                                  const char *default_instances, const int known_optima);
extern const size_t suite_biobj_ext_instances[][3];
