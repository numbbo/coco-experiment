#pragma once
#include "coco.h"
#include "coco_internal.h"
#undef INTERFACE
void observer_bbob_old(coco_observer_t *observer, const char *options, coco_option_keys_t **option_keys);
void logger_bbob_old_free(void *logger);
coco_problem_t *logger_bbob_old(coco_observer_t *observer, coco_problem_t *problem);
