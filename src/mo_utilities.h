#pragma once
#include "coco.h"
#undef INTERFACE
double mo_get_distance_to_ROI(const double *normalized_y, const size_t num_obj);
int mo_is_within_ROI(const double *normalized_y, const size_t num_obj);
int mo_get_dominance(const double *normalized_y1, const double *normalized_y2, const size_t num_obj);
double *mo_normalize(const double *y, const double *ideal, const double *nadir, const size_t num_obj);
double mo_get_norm(const double *first, const double *second, const size_t dim);
extern const double mo_discretization;
extern const double mo_precision;
