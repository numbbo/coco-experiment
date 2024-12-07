#pragma once
#include "coco.h"
#include "coco_internal.h"

#include <stdarg.h>

char *coco_string_trim(char *string);
size_t *coco_string_parse_ranges(const char *string, const size_t min, const size_t max, const char *name,
                                 const size_t max_count);
char *coco_remove_from_string(const char *string, const char *from, const char *to);
char **coco_string_split(const char *string, const char delimiter);
long coco_strfind(const char *base, const char *seq);
char *coco_strconcat(const char *s1, const char *s2);
char *coco_strdupf(const char *str, ...);
char *coco_vstrdupf(const char *str, va_list args);
char *coco_strdup(const char *string);
char *coco_allocate_string(const size_t number_of_elements);
size_t *coco_allocate_vector_size_t(const size_t number_of_elements);
