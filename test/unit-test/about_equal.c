#include "about_equal.h"

#include <stdio.h>
#include <math.h>

int about_equal_value(const double a, const double b) {
  const double absolute_error = fabs(a - b);
  const double larger = fabs(a) > fabs(b) ? a : b;
  const double relative_error = fabs((a - b) / (fabs(larger) + 1e-23));

  if (a == b || absolute_error < 1e-13 || relative_error < 4e-6) {
    return 1;
  }
  fprintf(stderr, "ERROR: about_equal_value failed: %f != %f (abs_err=%f, rel_err=%f)\n", a, b, absolute_error, relative_error);
  return 0;
}

int about_equal_vector(const double *a, const double *b, const size_t dimension) {
  size_t i;

  for (i = 0; i < dimension; i++) {
    if (!about_equal_value(a[i], b[i]))
      return 0;
  }
  return 1;
}

int about_equal_2d(const double *a, const double b1, const double b2) {
  return about_equal_value(a[0], b1) && about_equal_value(a[1], b2);
}
