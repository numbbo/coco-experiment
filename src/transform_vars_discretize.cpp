/**
 * @file transform_vars_discretize.cpp
 *
 * @brief Implementation of transforming a continuous problem to a mixed-integer problem by making some
 * of its variables discrete. The integer variables are considered as bounded (any variable outside the
 * decision space is mapped to the closest boundary point), while the continuous ones are treated as
 * unbounded.
 *
 * @note The first problem->number_of_integer_variables are integer, while the rest are continuous.
 *
 * The discretization works as follows. Consider the case where the interval [l, u] of the inner problem
 * needs to be discretized to n integer values of the outer problem. First, [l, u] is discretized to n
 * integers by placing the integers so that there is a (u-l)/(n+1) distance between them (and the border
 * points). Then, the transformation is shifted so that the optimum aligns with the closest integer. In
 * this way, we make sure that if the optimum is within [l, u], so are all the shifted points.
 *
 * When evaluating such a problem, the x values of the integer variables are first discretized. Any value
 * x < 0 is mapped to 0 and any value x > (n-1) is mapped to (n-1).
 */

#include <assert.h>

#include "coco.h"
#include "coco_problem.cpp"

namespace coco {

/**
 * @brief Makes the first number_of_integer_variables variables of the inner problem discrete.
 */
class TransformVarsDiscretize final : public TransformedProblem {
public:
  TransformVarsDiscretize(ProblemPtr inner, std::span<const double> smallest_values_of_interest,
                          std::span<const double> largest_values_of_interest, const size_t number_of_integer_variables)
      : TransformedProblem(std::move(inner), "transform_vars_discretize"),
        offset_(inner_problem_->number_of_variables()), discretized_x_(inner_problem_->number_of_variables()) {
    double l, u, inner_l, inner_u, outer_l, outer_u;
    double outer_xopt, inner_xopt, inner_approx_xopt;
    const double precision_offset = 1e-7; /* Needed to avoid issues with rounding doubles */
    int n;

    assert(number_of_integer_variables > 0);
    assert(smallest_values_of_interest.size() == number_of_variables_ &&
           largest_values_of_interest.size() == number_of_variables_);
    this->number_of_integer_variables_ = number_of_integer_variables;

    for (size_t i = 0; i < number_of_variables_; i++) {
      assert(smallest_values_of_interest[i] < largest_values_of_interest[i]);
      this->smallest_values_of_interest_[i] = smallest_values_of_interest[i];
      this->largest_values_of_interest_[i] = largest_values_of_interest[i];
      offset_[i] = 0;
      if (i < number_of_integer_variables) {
        /* Compute the offset for integer variables */
        outer_l = this->smallest_values_of_interest_[i];
        outer_u = this->largest_values_of_interest_[i];
        l = inner_problem_->smallest_values_of_interest()[i];
        u = inner_problem_->largest_values_of_interest()[i];
        n = DoubleToInt(outer_u) - DoubleToInt(outer_l) + 1; /* number of integer values */
        assert(n > 1);
        inner_l = l + (u - l) / (n + 1);
        inner_u = u - (u - l) / (n + 1);
        /* Find the location of the optimum in the coordinates of the outer problem */
        inner_xopt = inner_problem_->best_parameter()[i];
        outer_xopt = outer_l + (outer_u - outer_l) * (inner_xopt - inner_l) / (inner_u - inner_l);
        outer_xopt = DoubleRound(outer_xopt + precision_offset);
        /* Make sure you the bounds are respected */
        if (outer_xopt < outer_l)
          outer_xopt = outer_l;
        if (outer_xopt > outer_u)
          outer_xopt = outer_u;
        best_parameter_[i] = outer_xopt;
        /* Find the value corresponding to outer_xopt in the coordinates of the inner problem */
        inner_approx_xopt = inner_l + (inner_u - inner_l) * (outer_xopt - outer_l) / (outer_u - outer_l);
        /* Compute the difference between the inner_approx_xopt and inner_xopt */
        offset_[i] = inner_approx_xopt - inner_xopt;
      }
    }

    if (number_of_constraints_ > 0)
      coco_error("TransformVarsDiscretize(): Constraints not supported yet.");
  }

  /**
   * @brief Evaluates the transformed objective function.
   */
  void EvaluateFunction(std::span<const double> x, std::span<double> y) override {
    double l, u, inner_l, inner_u, outer_l, outer_u;
    int n;

    if (inner_problem_->number_of_objectives() == 0) {
      inner_problem_->EvaluateFunction(x, y);
      return;
    }

    /* Transform x to fit in the discretized space */
    discretized_x_.assign(x.begin(), x.end());
    for (size_t i = 0; i < number_of_integer_variables_; ++i) {
      outer_l = smallest_values_of_interest_[i];
      outer_u = largest_values_of_interest_[i];
      l = inner_problem_->smallest_values_of_interest()[i];
      u = inner_problem_->largest_values_of_interest()[i];
      n = DoubleToInt(outer_u) - DoubleToInt(outer_l) + 1; /* number of integer values in this coordinate */
      assert(n > 1);
      inner_l = l + (u - l) / (n + 1);
      inner_u = u - (u - l) / (n + 1);
      /* Make sure you the bounds are respected */
      discretized_x_[i] = DoubleRound(x[i]);
      if (discretized_x_[i] < outer_l)
        discretized_x_[i] = outer_l;
      if (discretized_x_[i] > outer_u)
        discretized_x_[i] = outer_u;
      discretized_x_[i] =
          inner_l + (inner_u - inner_l) * (discretized_x_[i] - outer_l) / (outer_u - outer_l) - offset_[i];
    }

    inner_problem_->EvaluateFunction(discretized_x_, y);
  }

  /** @brief Constraints are not supported yet (TODO?). */
  void EvaluateConstraint(std::span<const double> x, std::span<double> y, int update_counter) override {
    Problem::EvaluateConstraint(x, y, update_counter);
  }

  /** @brief Gradients are not supported yet (TODO?). */
  void EvaluateGradient(std::span<const double> x, std::span<double> y) override { Problem::EvaluateGradient(x, y); }

  /** @brief Rounding to the integers is not invertible. */
  bool Inverse(std::span<const double>, std::span<double>) const override { return false; }

private:
  std::vector<double> offset_;
  std::vector<double> discretized_x_;
};

} // namespace coco
