#pragma once
#include "coco.h"

typedef struct f_griewank_rosenbrock_data_s f_griewank_rosenbrock_data_t;

coco_problem_t *f_griewank_rosenbrock_permblockdiag_bbob_bbob_problem_allocate(
    const size_t function, const size_t dimension, const size_t instance, const long rseed, const void *args,
    const char *problem_id_template, const char *problem_name_template);
coco_problem_t *f_griewank_rosenbrock_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                            const size_t instance, const long rseed, const void *args,
                                                            const char *problem_id_template,
                                                            const char *problem_name_template);
coco_problem_t *f_griewank_rosenbrock_allocate(const size_t number_of_variables, double facftrue);
void f_griewank_rosenbrock_evaluate(coco_problem_t *problem, const double *x, double *y);
double f_griewank_rosenbrock_raw(const double *x, const size_t number_of_variables, f_griewank_rosenbrock_data_t *data);
