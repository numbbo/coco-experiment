#pragma once
#include "coco.h"
#include "coco_internal.h"
#undef INTERFACE
double coco_random_normal(coco_random_state_t *state);
double coco_random_uniform(coco_random_state_t *state);
void coco_random_free(coco_random_state_t *state);
coco_random_state_t *coco_random_new(uint32_t seed);
void coco_random_generate(coco_random_state_t *state);
