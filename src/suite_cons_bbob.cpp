/**
 * @file  suite_cons_bbob.cpp
 * @brief Implementation of the constrained bbob suite containing
 *        54 constrained problems in 6 dimensions. See comments in
 *        "suite_cons_bbob_problems.cpp" for more details.
 */

#include "coco.h"
#include "suite_cons_bbob_problems.cpp"
#include "transform_obj_scale.cpp"

namespace coco {

/**
 * @brief Sets the instances associated with years for the constrained
 *        bbob suite.
 */
static std::string_view SuiteConsBbobGetInstancesByYear(const int year) {

  if ((year >= 2022) || (year == 0)) {
    return "1-15";
  } else if ((year >= 2009)) {
    return "1-15";
  } else {
    coco_error("SuiteConsBbobGetInstancesByYear(): year %d not defined for suite_cons_bbob", year);
    return {};
  }
}

/**
 * @brief Creates and returns a constrained BBOB problem.
 */
static ProblemPtr GetConsBbobProblem(std::string_view suite_name, const size_t function, const size_t dimension,
                                     const size_t instance) {

  ProblemPtr problem;

  std::vector<double> feasible_direction(dimension);
  std::vector<double> xopt(dimension);
  long rseed = static_cast<long>(function + 10000 * instance);

  const char *problem_id_template = "bbob-constrained_f%03zu_i%02zu_d%02zu";
  const char *problem_name_template = "bbob-constrained suite problem f%zu instance %zu in %zuD";
  if (suite_name == "bbob-constrained-active-only") {
    /* CAVEAT: the first 28 chars of this ID are used in CLinearConsBbobProblemAllocate()
     * in c_linear.cpp to make the decision whether inactive constraints are added */
    problem_id_template = "bbob-constrained-active-only_f%03zu_i%02zu_d%02zu";
    problem_name_template = "bbob-constrained-active-only suite problem f%zu instance %zu in %zuD";
  }
  if (suite_name == "bbob-constrained-no-disguise") {
    /* First trick as above, by chance same number of characters
     * Both options at the same time not supported */
    problem_id_template = "bbob-constrained-no-disguise_f%03zu_i%02zu_d%02zu";
    problem_name_template = "bbob-constrained-no-disguise suite problem f%zu instance %zu in %zuD";
  }

  /* Seed value used for shifting the whole constrained problem */
  Bbob2009ComputeXopt(xopt, rseed);

  /* Choose a different seed value for building the objective function */
  rseed = static_cast<long>(function + 20000 * instance);

  size_t number_of_linear_constraints = NbOfLinearConstraints(function, dimension);

  if (ObjFunctionType(function) == 1) {

    problem =
        SphereCLinearConsBbobProblemAllocate(function, dimension, instance, number_of_linear_constraints, rseed,
                                             feasible_direction, xopt, problem_id_template, problem_name_template);
    problem = std::make_unique<TransformObjScale>(
        std::move(problem), 10.); /* move initial feasible point to a delta-f between 100 and 10000 */

  } else if (ObjFunctionType(function) == 2) {

    problem =
        EllipsoidCLinearConsBbobProblemAllocate(function, dimension, instance, number_of_linear_constraints, rseed,
                                                feasible_direction, xopt, problem_id_template, problem_name_template);
    problem = std::make_unique<TransformObjScale>(
        std::move(problem), 1e-4); /* move initial feasible point to a delta-f between 100 and 10000 */

  } else if (ObjFunctionType(function) == 3) {

    problem =
        LinearSlopeCLinearConsBbobProblemAllocate(function, dimension, instance, number_of_linear_constraints, rseed,
                                                  feasible_direction, xopt, problem_id_template, problem_name_template);
    problem = std::make_unique<TransformObjScale>(
        std::move(problem), 10.); /* move initial feasible point to a delta-f between 100 and 10000 */

  } else if (ObjFunctionType(function) == 4) {

    problem = EllipsoidRotatedCLinearConsBbobProblemAllocate(function, dimension, instance,
                                                             number_of_linear_constraints, rseed, feasible_direction,
                                                             xopt, problem_id_template, problem_name_template);
    problem = std::make_unique<TransformObjScale>(
        std::move(problem), 1e-4); /* move initial feasible point to a delta-f between 100 and 10000 */

  } else if (ObjFunctionType(function) == 5) {

    problem =
        DiscusCLinearConsBbobProblemAllocate(function, dimension, instance, number_of_linear_constraints, rseed,
                                             feasible_direction, xopt, problem_id_template, problem_name_template);
    problem = std::make_unique<TransformObjScale>(
        std::move(problem), 1e-4); /* move initial feasible point to a delta-f between 100 and 10000 */

  } else if (ObjFunctionType(function) == 6) {

    problem =
        BentCigarCLinearConsBbobProblemAllocate(function, dimension, instance, number_of_linear_constraints, rseed,
                                                feasible_direction, xopt, problem_id_template, problem_name_template);
    problem = std::make_unique<TransformObjScale>(
        std::move(problem), 1e-4); /* move initial feasible point to a delta-f between 100 and 100,000 */

  } else if (ObjFunctionType(function) == 7) {

    problem = DifferentPowersCLinearConsBbobProblemAllocate(function, dimension, instance, number_of_linear_constraints,
                                                            rseed, feasible_direction, xopt, problem_id_template,
                                                            problem_name_template);
    problem = std::make_unique<TransformObjScale>(
        std::move(problem), 1e2); /* move initial feasible point to a delta-f between 100 and 1000 */

  } else if (ObjFunctionType(function) == 8) {

    problem =
        RastriginCLinearConsBbobProblemAllocate(function, dimension, instance, number_of_linear_constraints, rseed,
                                                feasible_direction, xopt, problem_id_template, problem_name_template);
    problem = std::make_unique<TransformObjScale>(
        std::move(problem), 10.); /* move initial feasible point to a delta-f between 100 and 10000 */

  } else if (ObjFunctionType(function) == 9) {

    problem = RastriginRotatedCLinearConsBbobProblemAllocate(function, dimension, instance,
                                                             number_of_linear_constraints, rseed, feasible_direction,
                                                             xopt, problem_id_template, problem_name_template);
    problem = std::make_unique<TransformObjScale>(
        std::move(problem), 10.); /* move initial feasible point to a delta-f between 100 and 10000 */

  } else {
    coco_error("GetConsBbobProblem(): cannot retrieve problem f%zu instance %zu in %zuD", function, instance,
               dimension);
    return nullptr; /* Never reached */
  }

  return problem;
}

/**
 * @brief The bbob-constrained suites.
 */
class SuiteConsBbob final : public Suite {
public:
  using Suite::Suite;

  std::string_view GetInstancesByYear(const int year) const override { return SuiteConsBbobGetInstancesByYear(year); }

  /** @brief Returns the problem with the given function, dimension and instance. */
  ProblemPtr GetProblem(const size_t function, const size_t dimension, const size_t instance) override {
    ProblemPtr problem;

    problem = GetConsBbobProblem(name(), function, dimension, instance);

    /* Use the standard stacked problem_id as problem_name and
     * construct a new suite-specific problem_id
     */
    problem->SetId("%s_f%03zu_i%02zu_d%02zu", name().c_str(), function, instance, dimension);
    problem->SetName("%s suite problem f%zu instance %zu in %zuD", name().c_str(), function, instance, dimension);

    return problem;
  }
};

/**
 * @brief Sets the dimensions and default instances for the bbob suite.
 */
static Suite *SuiteConsBbobInitialize(std::string_view suite_name) {

  Suite *suite;
  const size_t dimensions[] = {2, 3, 5, 10, 20, 40};
  const size_t num_dimensions = sizeof(dimensions) / sizeof(dimensions[0]);

  /* IMPORTANT: Make sure to change the default instance for every new workshop! */
  suite = new SuiteConsBbob(suite_name, 54, num_dimensions, dimensions, "year: 2022", 1);

  return suite;
}

} // namespace coco
