#pragma once
#include "coco.h"

/**
 * @brief Data type for the versatile_data_t
 */
typedef struct {
  size_t proportion_short_axes_denom;
} f_discus_generalized_versatile_data_t;

coco_problem_t *f_discus_generalized_permblockdiag_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                         const size_t instance, const long rseed,
                                                                         const char *problem_id_template,
                                                                         const char *problem_name_template);
coco_problem_t *f_discus_generalized_allocate(const size_t number_of_variables, size_t proportion_short_axes_denom);
void f_discus_generalized_evaluate(coco_problem_t *problem, const double *x, double *y);
double f_discus_generalized_raw(const double *x, const size_t number_of_variables,
                                f_discus_generalized_versatile_data_t *f_discus_generalized_versatile_data);
void f_discus_generalized_versatile_data_free(coco_problem_t *problem);
