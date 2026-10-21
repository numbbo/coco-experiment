/**
 * @file f_schwefel_generalized.cpp
 * @brief Implementation of the Schwefel function of the large-scale suite.
 */

#include <math.h>

#include "coco.h"
#include "coco_internal.h"

namespace coco {

/**
 * @brief Implements the Schwefel function for large scale.
 */
static double SchwefelGeneralizedRaw(std::span<const double> x) {
  const size_t number_of_variables = x.size();

  size_t i = 0;
  double penalty, sum;
  const double Schwefel_constant = 4.189828872724339;

  /* Boundary handling*/
  penalty = 0.0;
  for (i = 0; i < number_of_variables; ++i) {
    const double tmp = fabs(x[i]) - 500.0;
    if (tmp > 0.0)
      penalty += tmp * tmp;
  }

  /* Computation core */
  sum = 0.0;
  for (i = 0; i < number_of_variables; ++i) {
    sum += x[i] * sin(sqrt(fabs(x[i])));
  }
  double result = 0.01 * (penalty + Schwefel_constant * 100. - sum / static_cast<double>(number_of_variables));

  return result;
}

} // namespace coco
