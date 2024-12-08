#pragma once
#include "coco.h"
#include "coco_internal.h"

extern coco_log_level_type_e coco_log_level;

int coco_vector_is_zero(const double *x, const size_t dim);
int coco_is_orthogonal(const double *M, const size_t nb_rows, const size_t nb_columns);
double coco_vector_scalar_product(const double *x, const double *y, size_t dimension);
double coco_vector_norm(const double *x, size_t dimension);
double coco_vector_scale(double *x, size_t dimension, double nom, double denom);
size_t coco_count_numbers(const size_t *numbers, const size_t max_count, const char *name);
char *coco_current_time_get_string(void);
int coco_is_feasible(coco_problem_t *problem, const double *x, double *constraint_values);
int coco_vector_isfinite(const double *x, const size_t dim);
void coco_vector_set_to_inf(double *y, const size_t dim);
int coco_is_inf(const double x);
void coco_vector_set_to_nan(double *y, const size_t dim);
int coco_vector_contains_nan(const double *x, const size_t dim);
int coco_is_nan(const double x);
int coco_double_to_int(const double number);
size_t coco_double_to_size_t(const double number);
double coco_double_min(const double a, const double b);
double coco_double_max(const double a, const double b);
double coco_double_round_up_with_precision(const double number, const double precision);
double coco_double_round(const double number);
int coco_double_almost_equal(const double a, const double b, const double precision);
int coco_options_read_values(const char *options, const char *name, char *pointer);
int coco_options_read_string(const char *options, const char *name, char *pointer);
int coco_options_read_double(const char *options, const char *name, double *pointer);
int coco_options_read_size_t(const char *options, const char *name, size_t *pointer);
int coco_options_read_int(const char *options, const char *name, int *pointer);
int coco_options_read(const char *options, const char *name, const char *format, void *pointer);
char *coco_option_keys_get_output_string(const coco_option_keys_t *option_keys, const char *info_string);
coco_option_keys_t *coco_option_keys(const char *option_string);
void coco_option_keys_add(coco_option_keys_t **basic_option_keys, const coco_option_keys_t *additional_option_keys);
coco_option_keys_t *coco_option_keys_get_redundant(const coco_option_keys_t *known_option_keys,
                                                   const coco_option_keys_t *given_option_keys);
void coco_option_keys_free(coco_option_keys_t *option_keys);
coco_option_keys_t *coco_option_keys_allocate(const size_t count, const char **keys);
double *coco_duplicate_vector(const double *src, const size_t number_of_elements);
char *coco_allocate_string(const size_t number_of_elements);
size_t *coco_allocate_vector_size_t(const size_t number_of_elements);
double *coco_allocate_vector_with_value(const size_t number_of_elements, double value);
double *coco_allocate_vector(const size_t number_of_elements);
int coco_remove_file(const char *file_name);
int coco_remove_directory(const char *path);
void coco_create_unique_directory(char **path);
void coco_create_unique_filename(const char *path, char **file_name, const char *ext);
int coco_create_directory(const char *path);
int coco_mkdir(const char *path);
int coco_file_exists(const char *path);
int coco_directory_exists(const char *path);
void coco_join_path(char *path, const size_t path_max_length, ...);
const char *coco_set_log_level(const char *log_level);
