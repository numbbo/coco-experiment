#include "minunit.h"

#include "coco.cpp"

namespace coco {

static double my_sqrt(double x) {
  return sqrt(x);
}

MU_TEST(test_brentq) {
  double y;

  y = Brentq(my_sqrt, 4, 1, 4, 1E-14, 1E-10, 200);
  mu_assert_double_eq(2.0, y);
}

MU_TEST(test_brentinv) {
  double y;

  y = BrentInv([](double x) { return x * x * x; }, -8);
  mu_assert_double_eq(-2.0, y);
}

MU_TEST(test_brentinv_function) {
  double y;

  y = BrentInv(my_sqrt, 4);
  mu_assert_double_eq(16.0, y);
}

} // namespace coco

/**
 * Run all tests in this file.
 */
int main(void) {
  MU_RUN_TEST(coco::test_brentq);
  MU_RUN_TEST(coco::test_brentinv);
  MU_RUN_TEST(coco::test_brentinv_function);

  MU_REPORT();

  return MU_EXIT_CODE;
}
