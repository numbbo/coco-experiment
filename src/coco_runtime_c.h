#pragma once
#include "coco.h"
#include "coco_internal.h"
#undef INTERFACE
void coco_free_memory(void *data);
void *coco_allocate_memory(const size_t size);
void coco_debug(const char *message, ...);
void coco_info_partial(const char *message, ...);
void coco_info(const char *message, ...);
void coco_warning(const char *message, ...);
void coco_error(const char *message, ...);
