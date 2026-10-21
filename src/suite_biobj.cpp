/**
 * @file suite_biobj.cpp
 * @brief Implementation of two bi-objective suites created by combining two single-objective problems
 * from the bbob suite:
 * - bbob-biobj contains 55 functions and 6 dimensions
 * - bbob-biobj-ext contains 55 + 37 functions and 6 dimensions
 *
 * The 55 functions of the bbob-biobj suite are created by combining any two single-objective bbob functions
 * i,j (where i<j) from a subset of 10 functions.
 *
 * The first 55 functions of the bbob-biobj-ext suite are the same as in the original bbob-biobj test suite
 * to which 37 functions are added. Those additional functions are constructed by combining all not yet
 * contained in-group combinations (i,j) of single-objective bbob functions i and j such that i<j (i.e. in
 * particular not all combinations (i,i) are included in this bbob-biobj-ext suite), with the exception of
 * the Weierstrass function (f16) for which the optimum is not unique and thus a nadir point is difficult
 * to compute, see http://numbbo.github.io/coco-doc/bbob-biobj/functions/ for details.
 *
 * @note See file suite_biobj_utilities.cpp for the implementation of the bi-objective problems and the handling
 * of new instances.
 */

#include "coco.h"
#include "mo_utilities.cpp"
#include "suite_biobj_utilities.cpp"
#include "suite_bbob.cpp"

namespace coco {

/**
 * @brief The bbob-biobj and bbob-biobj-ext suites.
 */
class SuiteBiobj final : public Suite {
public:
  /** @brief Creates the bbob-biobj or (with more functions) the bbob-biobj-ext suite. */
  explicit SuiteBiobj(std::string_view suite_name)
      : Suite(suite_name, suite_name == "bbob-biobj-ext" ? 55 + 37 : 55, {2, 3, 5, 10, 20, 40}, "instances: 1-15",
              true) {}

  /** @brief The instances of the bi-objective suites generally do not change with years. */
  std::string_view GetInstancesByYear(const int year) const override {
    if ((year == 2016) || (year == 0000)) { /* test case */
      return "1-10";
    } else
      return "1-15";
  }

  /** @brief Returns the problem with the given function, dimension and instance. */
  ProblemPtr GetProblem(const size_t function, const size_t dimension, const size_t instance) override {
    ProblemPtr problem;

    problem = GetBiobjProblem(function, dimension, instance, GetBbobProblem, new_instances_, dimensions().data(),
                              dimensions().size());

    problem->set_opt_known(true);

    return problem;
  }

private:
  SuiteBiobjNewInstances new_instances_; /**< @brief The new instances computed so far. */
};

} // namespace coco
