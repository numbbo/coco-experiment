#pragma once
#include "coco.h"
#undef INTERFACE
double coco_boundary_handling(coco_problem_t *problem, const double *x);
double coco_sample_uniform_noise(void);
double coco_sample_gaussian_noise(void);
void coco_reset_seeds(void);
void increase_random_seed(void);
void increase_random_n_seed(void);
extern long RANDNSEED;
extern long RANDSEED;
