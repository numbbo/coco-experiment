/**
 * @file bbob_problem.cpp
 * @brief The base class of the problems of the bbob, bbob-boxed and bbob-largescale suites.
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"
#include "suite_bbob_legacy_code.cpp"
#include "transform_functions.cpp"

/**
 * @brief The base class of the single-objective bbob problems, which apply their transformations explicitly.
 *
 * Derived implements double operator()(std::span<const double> x), the objective value at x, and sets the best
 * parameter and the best value in its constructor. The region of interest is [-5, 5]^n.
 */
template <typename Derived> class bbob_problem_t : public coco_problem_s {
public:
  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    y[0] = static_cast<Derived &>(*this)(x);
    assert(y[0] + 1e-13 >= best_value[0]);
  }

protected:
  /**
   * @brief Sets the region of interest and the problem id, name and type; the templates are formatted with the
   * function, instance and dimension.
   */
  bbob_problem_t(const size_t function, const size_t dimension, const size_t instance, const char *type,
                 const char *problem_id_template, const char *problem_name_template)
      : coco_problem_s(dimension, 1, 0) {
    init_from_scalars("", -5.0, 5.0, 0.0);
    set_id(problem_id_template, function, instance, dimension);
    set_name(problem_name_template, function, instance, dimension);
    set_type("%s", type);
  }

  /** @brief Sets the best value to the objective value at the best parameter. */
  void set_best_value_at_best_parameter() { best_value[0] = static_cast<Derived &>(*this)(best_parameter); }
};
