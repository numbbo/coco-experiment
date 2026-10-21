/**
 * @file bbob_problem.cpp
 * @brief The base class of the problems of the bbob, bbob-boxed and bbob-largescale suites.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_functions.cpp"

namespace coco {

/**
 * @brief The base class of the single-objective bbob problems, which apply their transformations explicitly.
 *
 * Derived implements double operator()(std::span<const double> x), the objective value at x, and sets the best
 * parameter and the best value in its constructor. The region of interest is [-5, 5]^n.
 */
template <typename Derived> class BbobProblem : public Problem {
public:
  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    y[0] = static_cast<Derived &>(*this)(x);
    assert(y[0] + 1e-13 >= best_value_[0]);
  }

protected:
  /**
   * @brief Sets the region of interest and the problem id, name and type; the templates are formatted with the
   * function, instance and dimension.
   */
  BbobProblem(const size_t function, const size_t dimension, const size_t instance, const char *type,
              const char *problem_id_template, const char *problem_name_template)
      : Problem(dimension, 1, 0) {
    InitFromScalars("", -5.0, 5.0, 0.0);
    SetId(problem_id_template, function, instance, dimension);
    SetName(problem_name_template, function, instance, dimension);
    SetType("%s", type);
  }

  /** @brief Sets the best value to the objective value at the best parameter. */
  void SetBestValueAtBestParameter() { best_value_[0] = static_cast<Derived &>(*this)(best_parameter_); }
};

} // namespace coco
