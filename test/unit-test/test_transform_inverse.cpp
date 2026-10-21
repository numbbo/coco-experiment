#include "minunit.h"

#include "coco.cpp"

static const size_t dimension = 6;

/**
 * A base problem that remembers the decision values it was last evaluated at.
 */
struct recorder_t final : public coco_problem_s {
  std::vector<double> last_x;

  recorder_t(const size_t number_of_objectives, const size_t number_of_constraints)
      : coco_problem_s(dimension, number_of_objectives, number_of_constraints), last_x(dimension) {
    for (size_t i = 0; i < dimension; ++i) {
      smallest_values_of_interest[i] = -5.0;
      largest_values_of_interest[i] = 5.0;
      best_parameter[i] = 0.0;
    }
    if (number_of_objectives > 0)
      best_value[0] = 0.0;
    set_id("recorder_d%02lu", (unsigned long)dimension);
  }

  void evaluate_function(std::span<const double> x, std::span<double> y) override {
    last_x.assign(x.begin(), x.end());
    y[0] = 0.0;
  }

  void evaluate_constraint(std::span<const double> x, std::span<double> y, int update_counter) override {
    (void)update_counter;
    last_x.assign(x.begin(), x.end());
    y[0] = 0.0;
  }
};

/** Some point in the coordinates of the base problem, with both signs and different magnitudes. */
static std::vector<double> base_point() {
  std::vector<double> base_x(dimension);
  for (size_t i = 0; i < dimension; ++i)
    base_x[i] = 3.0 * sin(1.0 + (double)i) * (i % 3 == 0 ? 0.1 : 1.0);
  return base_x;
}

/**
 * Checks that evaluating the problem built by make at from_base_coordinates(base_x) passes base_x to the base
 * problem again, up to the given relative precision (0 means exactly).
 */
template <typename Make> static int maps_back(Make make, const double precision, const bool constraint = false) {
  auto owned = std::make_unique<recorder_t>(constraint ? 0 : 1, constraint ? 1 : 0);
  recorder_t *recorder = owned.get();
  coco_problem_ptr_t problem = make(std::move(owned));
  const std::vector<double> base_x = base_point();
  std::vector<double> x(dimension), y(1);

  if (!problem->from_base_coordinates(base_x, x))
    return 0;
  if (constraint)
    problem->evaluate_constraint(x, y, 0);
  else
    problem->evaluate_function(x, y);
  for (size_t i = 0; i < dimension; ++i) {
    if (fabs(recorder->last_x[i] - base_x[i]) > precision * coco_double_max(1.0, fabs(base_x[i])))
      return 0;
  }
  return 1;
}

/** Checks that from_base_coordinates() reports that the problem built by make cannot be inverted. */
template <typename Make> static int is_not_invertible(Make make) {
  coco_problem_ptr_t problem = make(std::make_unique<recorder_t>(1, 0));
  const std::vector<double> base_x = base_point();
  std::vector<double> x(dimension);
  return !problem->from_base_coordinates(base_x, x);
}

/**
 * Checks that the inverse of the transformation maps the image of some point back to the point, up to the given
 * relative precision (0 means exactly).
 */
template <typename Transformation> static int inverts(const Transformation &transformation, const double precision) {
  const std::vector<double> x = base_point();
  std::vector<double> y(dimension), inverse_y(dimension);
  transformation.apply(x, y);
  transformation.inverse(y, inverse_y);
  for (size_t i = 0; i < dimension; ++i) {
    if (fabs(inverse_y[i] - x[i]) > precision * coco_double_max(1.0, fabs(x[i])))
      return 0;
  }
  return 1;
}

MU_TEST(test_exact_inverses) {
  const std::vector<double> offset = {0.5, -1.0, 2.0, -0.25, 3.0, 0.0};
  const size_t permutation[] = {3, 0, 5, 1, 4, 2};

  mu_check(maps_back([&](coco_problem_ptr_t p) { return transform_vars_shift(std::move(p), offset, 0); }, 1e-15));
  mu_check(maps_back([](coco_problem_ptr_t p) { return transform_vars_enforce_box(std::move(p)); }, 0));
  mu_check(maps_back([](coco_problem_ptr_t p) { return transform_obj_shift(std::move(p), 10.0); }, 0));

  mu_check(inverts(vars_scale_t(2.5, dimension), 1e-15));
  mu_check(inverts(vars_conditioning_t(100.0, dimension), 1e-15));
  mu_check(inverts(vars_permutation_t(permutation, false), 0));
  mu_check(inverts(vars_permutation_t(permutation, true), 0));
  mu_check(inverts(vars_x_hat_t(4711, dimension), 0));
  mu_check(inverts(vars_brs_t(dimension), 1e-15));
  mu_check(inverts(vars_z_hat_t(offset), 1e-15));
}

MU_TEST(test_rotation_inverses) {
  size_t number_of_blocks;
  const std::vector<size_t> block_sizes = coco_get_block_sizes(&number_of_blocks, dimension, "bbob-largescale");
  mu_check(inverts(vars_blockrotation_t(coco_compute_blockrotation(12345, block_sizes)), 1e-13));
  mu_check(inverts(vars_permblockdiag_t(dimension, 1, 2, 3), 1e-13));

  const coco_matrix_t R = bbob2009_compute_rotation(12345, dimension);
  const std::vector<double> b = {0.5, -1.0, 2.0, -0.25, 3.0, 0.0};
  mu_check(maps_back([&](coco_problem_ptr_t p) { return transform_vars_affine(std::move(p), R, b); }, 1e-13));

  /* A matrix that is not orthogonal is inverted by solving the linear system */
  coco_matrix_t M = R;
  M(0, 0) *= 2.0;
  mu_check(maps_back([&](coco_problem_ptr_t p) { return transform_vars_affine(std::move(p), M, b); }, 1e-12));

  /* A singular matrix is not invertible */
  for (size_t j = 0; j < dimension; ++j)
    M(1, j) = M(0, j);
  mu_check(is_not_invertible([&](coco_problem_ptr_t p) { return transform_vars_affine(std::move(p), M, b); }));
}

MU_TEST(test_numerical_inverses) {
  mu_check(maps_back([](coco_problem_ptr_t p) { return transform_vars_oscillate(std::move(p)); }, 1e-9));
  mu_check(maps_back([](coco_problem_ptr_t p) { return transform_vars_asymmetric(std::move(p), 0.2); }, 1e-9));
  /* The composition used by the constrained suite */
  mu_check(maps_back(
      [](coco_problem_ptr_t p) {
        return transform_vars_oscillate(transform_vars_asymmetric(std::move(p), 0.2));
      },
      1e-9));
}

MU_TEST(test_chain_inverse) {
  const std::vector<double> offset = {0.5, -1.0, 2.0, -0.25, 3.0, 0.0};
  const coco_matrix_t R = bbob2009_compute_rotation(12345, dimension);
  mu_check(maps_back(
      [&](coco_problem_ptr_t p) {
        p = transform_vars_affine(std::move(p), R, offset);
        p = transform_obj_shift(std::move(p), 1.0);
        p = transform_vars_asymmetric(std::move(p), 0.2);
        return transform_vars_shift(std::move(p), offset, 0);
      },
      1e-9));
}

MU_TEST(test_constraint_inverse) {
  const std::vector<double> gradient = {1.0, -2.0, 0.5, 3.0, -0.25, 1.5};
  mu_check(
      maps_back([&](coco_problem_ptr_t p) { return std::make_unique<c_linear_t>(std::move(p), gradient, 0.5, 2.0); },
                1e-14, true));

  const std::vector<double> zero_gradient = {1.0, 0.0, 0.5, 3.0, -0.25, 1.5};
  mu_check(is_not_invertible(
      [&](coco_problem_ptr_t p) { return std::make_unique<c_linear_t>(std::move(p), zero_gradient, 0.5, 2.0); }));
}

MU_TEST(test_rounding_is_not_invertible) {
  const std::vector<double> lower(dimension, -5.0), upper(dimension, 5.0);
  mu_check(is_not_invertible(
      [&](coco_problem_ptr_t p) { return transform_vars_discretize(std::move(p), lower, upper, dimension); }));
}

/**
 * Run all tests in this file.
 */
int main(void) {
  MU_RUN_TEST(test_exact_inverses);
  MU_RUN_TEST(test_rotation_inverses);
  MU_RUN_TEST(test_numerical_inverses);
  MU_RUN_TEST(test_chain_inverse);
  MU_RUN_TEST(test_constraint_inverse);
  MU_RUN_TEST(test_rounding_is_not_invertible);

  MU_REPORT();

  return MU_EXIT_CODE;
}
