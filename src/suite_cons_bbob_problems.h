#pragma once
#include "coco.h"
#include "coco_internal.h"
#undef INTERFACE
coco_problem_t *f_rastrigin_rotated_c_linear_cons_bbob_problem_allocate(
    const size_t function, const size_t dimension, const size_t instance, const size_t number_of_linear_constraints,
    const long rseed, double *feasible_direction, const double *xopt, const char *problem_id_template,
    const char *problem_name_template);
coco_problem_t *f_rastrigin_c_linear_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                const size_t instance,
                                                                const size_t number_of_linear_constraints,
                                                                const long rseed, double *feasible_direction,
                                                                const double *xopt, const char *problem_id_template,
                                                                const char *problem_name_template);
coco_problem_t *f_different_powers_c_linear_cons_bbob_problem_allocate(
    const size_t function, const size_t dimension, const size_t instance, const size_t number_of_linear_constraints,
    const long rseed, double *feasible_direction, const double *xopt, const char *problem_id_template,
    const char *problem_name_template);
coco_problem_t *f_bent_cigar_c_linear_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                 const size_t instance,
                                                                 const size_t number_of_linear_constraints,
                                                                 const long rseed, double *feasible_direction,
                                                                 const double *xopt, const char *problem_id_template,
                                                                 const char *problem_name_template);
coco_problem_t *f_discus_c_linear_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                             const size_t instance,
                                                             const size_t number_of_linear_constraints,
                                                             const long rseed, double *feasible_direction,
                                                             const double *xopt, const char *problem_id_template,
                                                             const char *problem_name_template);
coco_problem_t *f_linear_slope_c_linear_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                   const size_t instance,
                                                                   const size_t number_of_linear_constraints,
                                                                   const long rseed, double *feasible_direction,
                                                                   const double *xopt, const char *problem_id_template,
                                                                   const char *problem_name_template);
coco_problem_t *f_ellipsoid_rotated_c_linear_cons_bbob_problem_allocate(
    const size_t function, const size_t dimension, const size_t instance, const size_t number_of_linear_constraints,
    const long rseed, double *feasible_direction, const double *xopt, const char *problem_id_template,
    const char *problem_name_template);
coco_problem_t *f_ellipsoid_c_linear_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                                const size_t instance,
                                                                const size_t number_of_linear_constraints,
                                                                const long rseed, double *feasible_direction,
                                                                const double *xopt, const char *problem_id_template,
                                                                const char *problem_name_template);
coco_problem_t *f_sphere_c_linear_cons_bbob_problem_allocate(const size_t function, const size_t dimension,
                                                             const size_t instance,
                                                             const size_t number_of_linear_constraints,
                                                             const long rseed, double *feasible_direction,
                                                             const double *xopt, const char *problem_id_template,
                                                             const char *problem_name_template);
void feasible_direction_set_length(double *feasible_direction, const double *xopt, size_t dimension, long rseed);
size_t nb_of_linear_constraints(const size_t function, const size_t dimension);
size_t obj_function_type(const size_t function);
