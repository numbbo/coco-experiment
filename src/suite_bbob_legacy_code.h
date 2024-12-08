#pragma once
#include "coco.h"

double bbob2009_compute_fopt(const size_t function, const size_t instance);
void sbox_cost_compute_xopt(double *xopt, const long seed, const size_t DIM);
void bbob2009_compute_xopt(double *xopt, const long seed, const size_t DIM);
void bbob2009_copy_rotation_matrix(double **rot, double *M, double *b, const size_t DIM);
void bbob2009_compute_rotation(double **B, const long seed, const size_t DIM);
void bbob2009_gauss(double *g, const size_t N, const long seed);
double **bbob2009_reshape(double **B, double *vector, const size_t m, const size_t n);
void bbob2009_unif(double *r, size_t N, long inseed);
void bbob2009_free_matrix(double **matrix, const size_t n);
double **bbob2009_allocate_matrix(const size_t n, const size_t m);
double bbob2009_round(double x);
double bbob2009_fmax(double a, double b);
double bbob2009_fmin(double a, double b);
