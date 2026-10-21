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
 * @brief Sets the instances associated with years for the bbob-biobj-mixint suites.
 *
 * @note The instances of the bi-objective suites generally do not changes with years.
 */
static std::string_view SuiteBiobjMixintGetInstancesByYear(const int) { return "1-15"; }

/**
 * @brief Creates and returns a mixed-integer bi-objective bbob problem without needing the actual
 * bbob-mixint suite.
 *
 * The problem is constructed by first finding the underlying single-objective continuous problems,
 * then discretizing the problems, then scaling them to adjust their difficulty and finally stacking
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

  ProblemPtr problem_cont;

  ProblemPtr problem;
  ProblemPtr problem1;
  ProblemPtr problem2;
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
  problem_cont = GetBiobjProblem(function, dimension, instance, coco_get_problem_function, new_instances, dimensions,
                                 num_dimensions);
  assert(problem_cont != nullptr);
  /* Take the two single-objective problems out of the stacked problem and free only the latter */
  problem1 = ProblemAsStacked(problem_cont.get())->ReleaseProblem1();
  problem2 = ProblemAsStacked(problem_cont.get())->ReleaseProblem2();
  assert(problem1);
  assert(problem2);
  problem_cont.reset();

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

  const size_t function1 = ProblemGetSuiteDepFunction(problem1.get());
  const size_t function2 = ProblemGetSuiteDepFunction(problem2.get());
  const std::string problem_type1 = problem1->type();
  const std::string problem_type2 = problem2->type();

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
  problem = std::make_unique<StackedProblem>(std::move(problem1_mixint), std::move(problem2_mixint),
                                             smallest_values_of_interest, largest_values_of_interest);

  /* Use the standard stacked problem_id as problem_name and construct a new problem_id */
  problem->SetName("%s", problem->id().c_str());
  problem->SetId("bbob-biobj-mixint_f%03zu_i%02zu_d%03zu", function, instance, dimension);

  /* Construct problem type */
  problem->SetType("%s_%s", problem_type1.c_str(), problem_type2.c_str());

  return problem;
}
/**
 * @brief The bbob-biobj-mixint suite.
 */
class SuiteBiobjMixint final : public Suite {
public:
  using Suite::Suite;

  std::string_view GetInstancesByYear(const int year) const override {
    return SuiteBiobjMixintGetInstancesByYear(year);
  }

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

/**
 * @brief Sets the dimensions and default instances for the bbob-biobj-mixint suite.
 */
static Suite *SuiteBiobjMixintInitialize() {

  Suite *suite;
  const size_t dimensions[] = {5, 10, 20, 40, 80, 160};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  suite = new SuiteBiobjMixint("bbob-biobj-mixint", 92, num_dimensions, dimensions, "instances: 1-15", 0);

  return suite;
}

} // namespace coco
