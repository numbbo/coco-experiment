/**
 * @file suite_bbob_mixint.cpp
 * @brief Implementation of a suite with mixed-integer bbob problems. The functions are the same
 * as those in the bbob suite with 24 functions, but the large-scale implementations of the
 * functions are used instead of the original ones for dimensions over 40. Additionally, the
 * functions are scaled as given by suite_bbob_mixint_scaling_factors.
 */

#include "coco.h"
#include "suite_bbob.cpp"
#include "suite_largescale.cpp"
#include "transform_vars_discretize.cpp"
#include "transform_obj_scale.cpp"

namespace coco {

/**
 * @brief Factors used to scale separate functions in order to achieve similar difficulty.
 */
static const double SUITE_BBOB_MIXINT_SCALING_FACTORS[] = {
    1,    1e-3, 1e-1, 1e-1, 1,    /* f1 to f5 */
    1e-2, 1,    1e-2, 1e-2, 1e-3, /* f6 to f10 */
    1e-2, 1e-4, 1e-1, 1,    1e-1, /* f11 to f15 */
    1,    10,   1,    10,   1e-1, /* f16 to f20 */
    1,    1,    10,   1e-1        /* f21 to f24 */
};

/**
 * @brief Creates and returns a mixed-integer bbob problem without needing the actual bbob-mixint
 * suite.
 *
 * @param function Function
 * @param dimension Dimension
 * @param instance Instance
 * @param coco_get_problem_function The function that is used to access the continuous problem.
 * @return The problem that corresponds to the given parameters.
 */
static ProblemPtr GetBbobMixintProblem(const size_t function, const size_t dimension, const size_t instance,
                                       const GetProblemFunction coco_get_problem_function) {
  ProblemPtr problem;

  /* The cardinality of variables (0 = continuous variables should always come last) */
  const size_t variable_cardinality[] = {2, 4, 8, 16, 0};

  std::vector<double> smallest_values_of_interest(dimension);
  std::vector<double> largest_values_of_interest(dimension);
  std::string inner_problem_id;

  size_t j;
  size_t cardinality = 0;
  size_t num_integer = dimension;
  if (dimension % 5 != 0)
    coco_error("GetBbobMixintProblem(): dimension %zu not supported for suite_bbob_mixint", dimension);

  problem = coco_get_problem_function(function, dimension, instance);
  assert(problem != nullptr);

  /* Set the ROI of the outer problem according to the given cardinality of variables and the ROI of the
   * inner problem to [-4, 4] for variables that will be discretized */
  for (size_t i = 0; i < dimension; i++) {
    j = i / (dimension / 5);
    cardinality = variable_cardinality[j];
    if (cardinality == 0) {
      /* Continuous variables */
      /* Outer problem */
      smallest_values_of_interest[i] = -5;
      largest_values_of_interest[i] = 5;
      if (num_integer == dimension)
        num_integer = i;
    } else {
      /* Outer problem */
      smallest_values_of_interest[i] = 0;
      largest_values_of_interest[i] = static_cast<double>(cardinality) - 1;
      /* Inner problem */
      problem->SetRegionOfInterest(i, -4, 4);
    }
  }

  inner_problem_id = problem->id();

  problem = std::make_unique<TransformVarsDiscretize>(std::move(problem), smallest_values_of_interest,
                                                      largest_values_of_interest, num_integer);

  problem = std::make_unique<TransformObjScale>(std::move(problem), SUITE_BBOB_MIXINT_SCALING_FACTORS[function - 1]);

  problem->SetId("bbob-mixint_f%03zu_i%02zu_d%02zu", function, instance, dimension);
  problem->SetName("mixint(%s)", inner_problem_id.c_str());

  return problem;
}

/**
 * @brief The bbob-mixint suite.
 */
class SuiteBbobMixint final : public Suite {
public:
  SuiteBbobMixint() : Suite("bbob-mixint", 24, {5, 10, 20, 40, 80, 160}, "instances: 1-15", true) {}

  std::string_view GetInstancesByYear(const int) const override { return "1-15"; }

  /** @brief Returns the problem with the given function, dimension and instance. */
  ProblemPtr GetProblem(const size_t function, const size_t dimension, const size_t instance) override {
    ProblemPtr problem;
    const size_t dim_large_scale = 50; /* Switch to large-scale functions for dimensions over 50 */

    if (dimension < dim_large_scale)
      problem = GetBbobMixintProblem(function, dimension, instance, GetBbobProblem);
    else
      problem = GetBbobMixintProblem(function, dimension, instance, GetLargescaleProblem);

    return problem;
  }
};

} // namespace coco
