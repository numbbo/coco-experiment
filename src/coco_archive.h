#pragma once
#include "coco.h"

typedef struct coco_archive_s coco_archive_t;
typedef struct coco_archive_avl_item_s coco_archive_avl_item_t;

void coco_archive_free(coco_archive_t *archive);
double coco_archive_get_hypervolume(coco_archive_t *archive);
size_t coco_archive_get_number_of_solutions(coco_archive_t *archive);
const char *coco_archive_get_next_solution_text(coco_archive_t *archive);
void coco_archive_update(coco_archive_t *archive);
int coco_archive_add_solution(coco_archive_t *archive, const double y1, const double y2, const char *text);
coco_archive_t *coco_archive(const char *suite_name, const size_t function, const size_t dimension,
                             const size_t instance);
coco_archive_t *coco_archive_allocate(void);
int coco_archive_compare_by_last_objective(const coco_archive_avl_item_t *item1, const coco_archive_avl_item_t *item2,
                                           void *userdata);
void coco_archive_node_item_free(coco_archive_avl_item_t *item, void *userdata);
coco_archive_avl_item_t *coco_archive_node_item_create(const double *y, const double *ideal, const double *nadir,
                                                       const size_t num_obj, const char *text);
