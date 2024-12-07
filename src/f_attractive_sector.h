#pragma once
#include "coco.h"

typedef struct f_attractive_sector_data_s f_attractive_sector_data_t;

coco_problem_t *f_attractive_sector_permblockdiag_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                        const size_t instance, const long rseed,
                                                                        const char *problem_id_template,
                                                                        const char *problem_name_template);
coco_problem_t *f_attractive_sector_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                          const size_t instance, const long rseed,
                                                          const char *problem_id_template,
                                                          const char *problem_name_template);
coco_problem_t *f_attractive_sector_allocate(const size_t number_of_variables, const double *xopt);
void f_attractive_sector_free(coco_problem_t *problem);
void f_attractive_sector_evaluate(coco_problem_t *problem, const double *x, double *y);
double f_attractive_sector_raw(const double *x, const size_t number_of_variables, f_attractive_sector_data_t *data);
