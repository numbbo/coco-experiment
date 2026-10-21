#ifndef ABOUT_EQUAL_H
#define ABOUT_EQUAL_H
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

int about_equal_value(const double a, const double b);

int about_equal_vector(const double *a, const double *b, const size_t dimension);

int about_equal_2d(const double *a, const double b1, const double b2);

#ifdef __cplusplus
}
#endif

#endif
