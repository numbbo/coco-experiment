/**
 * @file suite_biobj_mixint.cpp
 * @brief Implementation of a bi-objective mixed-integer suite. The functions are the same as those
 * in the bbob-biobj-ext suite with 92 functions, but the large-scale implementations of the
 * functions are used instead of the original ones for dimensions over 40.
 */

#include "coco.h"
#include "mo_utilities.cpp"
#include "suite_biobj_utilities.cpp"
#include "suite_largescale.cpp"
#include "transform_vars_discretize.cpp"
#include "transform_obj_scale.cpp"
#include "suite_bbob_mixint.cpp"

namespace coco {

/**
 * @brief Creates and returns a mixed-integer bi-objective bbob problem without needing the actual
 * bbob-mixint suite.
 *
 * The problem is constructed by first finding the underlying single-objective continuous problems,
 * then discretizing the problems, then scaling them to adjust their difficulty and finally combining
 * them to get a bi-objective mixed-integer problem.
 *
 * @param function Function
 * @param dimension Dimension
 * @param instance Instance
 * @param coco_get_problem_function The function that is used to access the single-objective problem.
 * @param new_instances The new instances computed so far (new ones are added).
 * @param dimensions An array of dimensions to take into account when creating new instances.
 * @param num_dimensions The number of dimensions to take into account when creating new instances.
 * @return The problem that corresponds to the given parameters.
 */
static ProblemPtr GetBiobjMixintProblem(const size_t function, const size_t dimension, const size_t instance,
                                        const GetProblemFunction coco_get_problem_function,
                                        SuiteBiobjNewInstances &new_instances, const size_t *dimensions,
                                        const size_t num_dimensions) {

  ProblemPtr problem;
  ProblemPtr problem1_mixint;
  ProblemPtr problem2_mixint;

  std::vector<double> smallest_values_of_interest(dimension);
  std::vector<double> largest_values_of_interest(dimension);

  size_t j;
  size_t num_integer = dimension;
  /* The cardinality of variables (0 = continuous variables should always come last) */
  const size_t variable_cardinality[] = {2, 4, 8, 16, 0};

  if (dimension % 5 != 0)
    coco_error("GetBiobjMixintProblem(): dimension %zu not supported for suite_bbob_mixint", dimension);

  /* First, find the underlying single-objective continuous problems */
  auto [problem1, problem2] = GetBiobjObjectives(function, dimension, instance, coco_get_problem_function,
                                                 new_instances, dimensions, num_dimensions);

  /* Set the ROI of the outer problem according to the given cardinality of variables and the ROI of the
   * inner problems to [-4, 4] for variables that will be discretized */
  for (size_t i = 0; i < dimension; i++) {
    j = i / (dimension / 5);
    if (variable_cardinality[j] == 0) {
      /* Continuous variables */
      /* Outer problem */
      smallest_values_of_interest[i] = -100;
      largest_values_of_interest[i] = 100;
      if (num_integer == dimension)
        num_integer = i;
    } else {
      /* Outer problem */
      smallest_values_of_interest[i] = 0;
      largest_values_of_interest[i] = static_cast<double>(variable_cardinality[j]) - 1;
      /* Inner problems */
      problem1->SetRegionOfInterest(i, -4, 4);
      problem2->SetRegionOfInterest(i, -4, 4);
    }
  }

  const size_t function1 = problem1->suite_dep_function();
  const size_t function2 = problem2->suite_dep_function();

  /* Second, discretize the single-objective problems */
  problem1_mixint = std::make_unique<TransformVarsDiscretize>(std::move(problem1), smallest_values_of_interest,
                                                              largest_values_of_interest, num_integer);
  problem2_mixint = std::make_unique<TransformVarsDiscretize>(std::move(problem2), smallest_values_of_interest,
                                                              largest_values_of_interest, num_integer);

  /* Third, scale the objective values */
  problem1_mixint =
      std::make_unique<TransformObjScale>(std::move(problem1_mixint), SUITE_BBOB_MIXINT_SCALING_FACTORS[function1 - 1]);
  problem2_mixint =
      std::make_unique<TransformObjScale>(std::move(problem2_mixint), SUITE_BBOB_MIXINT_SCALING_FACTORS[function2 - 1]);

  /* Fourth, combine the problems in a bi-objective mixed-integer problem */
  problem = std::make_unique<BiobjProblem>(std::move(problem1_mixint), std::move(problem2_mixint),
                                           smallest_values_of_interest, largest_values_of_interest);

  /* The name is made of the ids of both problems */
  problem->SetId("bbob-biobj-mixint_f%03zu_i%02zu_d%03zu", function, instance, dimension);

  return problem;
}
/**
 * @brief The bbob-biobj-mixint suite.
 */
class SuiteBiobjMixint final : public Suite {
public:
  SuiteBiobjMixint() : Suite("bbob-biobj-mixint", 92, {5, 10, 20, 40, 80, 160}, "instances: 1-15", false) {}

  std::string_view GetInstancesByYear(const int) const override { return "1-15"; }

  /** @brief Returns the problem with the given function, dimension and instance. */
  ProblemPtr GetProblem(const size_t function, const size_t dimension, const size_t instance) override {
    ProblemPtr problem;
    const size_t dim_large_scale = 50; /* Switch to large-scale functions for dimensions over 50 */

    if (dimension < dim_large_scale)
      problem = GetBiobjMixintProblem(function, dimension, instance, GetBbobProblem, new_instances_,
                                      dimensions().data(), dimensions().size());
    else
      problem = GetBiobjMixintProblem(function, dimension, instance, GetLargescaleProblem, new_instances_,
                                      dimensions().data(), dimensions().size());

    return problem;
  }

private:
  SuiteBiobjNewInstances new_instances_; /**< @brief The new instances computed so far. */
};

} // namespace coco
