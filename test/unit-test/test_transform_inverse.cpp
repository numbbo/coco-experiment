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

  void evaluate_function(const double *x, double *y) override {
    last_x.assign(x, x + dimension);
    y[0] = 0.0;
  }

  void evaluate_constraint(const double *x, double *y, int update_counter) override {
    (void)update_counter;
    last_x.assign(x, x + dimension);
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

  if (!problem->from_base_coordinates(base_x.data(), x.data()))
    return 0;
  if (constraint)
    problem->evaluate_constraint(x.data(), y.data(), 0);
  else
    problem->evaluate_function(x.data(), y.data());
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
  return !problem->from_base_coordinates(base_x.data(), x.data());
}

MU_TEST(test_exact_inverses) {
  const std::vector<double> offset = {0.5, -1.0, 2.0, -0.25, 3.0, 0.0};
  const std::vector<double> signs = {1, -1, -1, 1, 1, -1};
  const size_t permutation[] = {3, 0, 5, 1, 4, 2};

  mu_check(maps_back([&](coco_problem_ptr_t p) { return transform_vars_shift(std::move(p), offset.data(), 0); },
                     1e-15));
  mu_check(maps_back([](coco_problem_ptr_t p) { return transform_vars_scale(std::move(p), 2.5); }, 1e-15));
  mu_check(maps_back([](coco_problem_ptr_t p) { return transform_vars_conditioning(std::move(p), 100.0); }, 1e-15));
  mu_check(maps_back(
      [&](coco_problem_ptr_t p) { return transform_vars_permutation(std::move(p), permutation, dimension); }, 0));
  mu_check(maps_back(
      [&](coco_problem_ptr_t p) { return transform_vars_inverse_permutation(std::move(p), permutation, dimension); },
      0));
  mu_check(maps_back([](coco_problem_ptr_t p) { return transform_vars_x_hat(std::move(p), 4711); }, 0));
  mu_check(maps_back([](coco_problem_ptr_t p) { return transform_vars_brs(std::move(p)); }, 1e-15));
  mu_check(maps_back([&](coco_problem_ptr_t p) { return transform_vars_z_hat(std::move(p), offset.data()); }, 1e-15));
  mu_check(maps_back([](coco_problem_ptr_t p) { return transform_vars_enforce_box(std::move(p)); }, 0));
  mu_check(maps_back([](coco_problem_ptr_t p) { return transform_obj_shift(std::move(p), 10.0); }, 0));

  f_lunacek_bi_rastrigin_shared_data_t lunacek_shared;
  lunacek_shared.x_hat.resize(dimension);
  mu_check(maps_back(
      [&](coco_problem_ptr_t p) { return transform_vars_x_hat_generic(std::move(p), signs.data(), &lunacek_shared); },
      0));
}

MU_TEST(test_rotation_inverses) {
  size_t number_of_blocks;
  const std::vector<size_t> block_sizes = coco_get_block_sizes(&number_of_blocks, dimension, "bbob-largescale");
  const coco_block_diagonal_matrix_t B = coco_compute_blockrotation(12345, block_sizes.data(), number_of_blocks);
  mu_check(maps_back([&](coco_problem_ptr_t p) { return transform_vars_blockrotation(std::move(p), B); }, 1e-13));

  const coco_matrix_t R = bbob2009_compute_rotation(12345, dimension);
  const std::vector<double> b = {0.5, -1.0, 2.0, -0.25, 3.0, 0.0};
  mu_check(maps_back([&](coco_problem_ptr_t p) { return transform_vars_affine(std::move(p), R, b); }, 1e-13));

  /* A matrix that is not orthogonal */
  coco_matrix_t M = R;
  M(0, 0) *= 2.0;
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
  const size_t permutation[] = {3, 0, 5, 1, 4, 2};
  mu_check(maps_back(
      [&](coco_problem_ptr_t p) {
        p = transform_vars_permutation(std::move(p), permutation, dimension);
        p = transform_vars_conditioning(std::move(p), 10.0);
        p = transform_obj_shift(std::move(p), 1.0);
        p = transform_vars_brs(std::move(p));
        return transform_vars_shift(std::move(p), offset.data(), 0);
      },
      1e-14));
}

MU_TEST(test_constraint_inverse) {
  const std::vector<double> gradient = {1.0, -2.0, 0.5, 3.0, -0.25, 1.5};
  mu_check(maps_back(
      [&](coco_problem_ptr_t p) { return std::make_unique<c_linear_t>(std::move(p), gradient.data(), 0.5, 2.0); },
      1e-14, true));

  const std::vector<double> zero_gradient = {1.0, 0.0, 0.5, 3.0, -0.25, 1.5};
  mu_check(is_not_invertible(
      [&](coco_problem_ptr_t p) { return std::make_unique<c_linear_t>(std::move(p), zero_gradient.data(), 0.5, 2.0); }));
}

MU_TEST(test_rounding_is_not_invertible) {
  f_step_ellipsoid_shared_data_t step_shared;
  mu_check(is_not_invertible(
      [&](coco_problem_ptr_t p) { return transform_vars_round_step(std::move(p), 10.0, &step_shared); }));

  const std::vector<double> lower(dimension, -5.0), upper(dimension, 5.0);
  mu_check(is_not_invertible([&](coco_problem_ptr_t p) {
    return transform_vars_discretize(std::move(p), lower.data(), upper.data(), dimension);
  }));
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

  return minunit_status;
}
