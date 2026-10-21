#include "minunit.h"

#include "coco.cpp"
#include "about_equal.h"

namespace coco {

/**
 * Tests the function coco_observer_targets_trigger.
 */
MU_TEST(test_coco_observer_targets_trigger) {

  ObserverTargets targets = ObserverTargets(0, 1e-3, 1, 1e-5);
  int update;

  update = targets.Trigger(99.99999);
  mu_check(update);
  mu_check(about_equal_value(targets.LastTarget(), 100));

  update = targets.Trigger(0.99999);
  mu_check(update);
  mu_check(about_equal_value(targets.LastTarget(), 1));

  update = targets.Trigger(0.899999);
  mu_check(update);
  mu_check(about_equal_value(targets.LastTarget(), 0.9));

  update = targets.Trigger(0.0199999);
  mu_check(update);
  mu_check(about_equal_value(targets.LastTarget(), 0.02));

  update = targets.Trigger(0.00000099999);
  mu_check(update);
  mu_check(about_equal_value(targets.LastTarget(), 1e-3));

  /* An improvement on the lin targets should not trigger an update if it is not an improvement overall */
  update = targets.Trigger(0.00999);
  mu_check(!update);

  update = targets.Trigger(0);
  mu_check(update);
  mu_check(about_equal_value(targets.LastTarget(), 0));

  update = targets.Trigger(-0.00100001);
  mu_check(update);
  mu_check(about_equal_value(targets.LastTarget(), -1.0e-3));

  update = targets.Trigger(-0.99999);
  mu_check(update);
  mu_check(about_equal_value(targets.LastTarget(), -0.999));

  update = targets.Trigger(-0.100099999);
  mu_check(!update);

}

/**
 * Tests the function coco_observer_log_targets_trigger.
 */
MU_TEST(test_coco_observer_log_targets_trigger) {

  ObserverLogTargets targets = ObserverLogTargets(10, 1e-8);
  int update;

  update = targets.Trigger(99.99999);
  mu_check(update);
  mu_check(about_equal_value(targets.value(), 100));

  update = targets.Trigger(10);
  mu_check(update);
  mu_check(about_equal_value(targets.value(), 10));

  update = targets.Trigger(1.00000);
  mu_check(update);
  mu_check(about_equal_value(targets.value(), 1));

  update = targets.Trigger(10);
  mu_check(!update);

  update = targets.Trigger(2);
  mu_check(!update);

  update = targets.Trigger(0.00001);
  mu_check(update);
  mu_check(about_equal_value(targets.value(), 0.00001));

  update = targets.Trigger(12);
  mu_check(!update);

  update = targets.Trigger(1e-8);
  mu_check(update);
  mu_check(about_equal_value(targets.value(), 1e-8));

  update = targets.Trigger(1e-9);
  mu_check(!update);

  update = targets.Trigger(-1e-8);
  mu_check(update);
  mu_check(about_equal_value(targets.value(), -1e-8));

  update = targets.Trigger(-1e-7);
  mu_check(update);
  mu_check(about_equal_value(targets.value(), -1e-7));

  update = targets.Trigger(2);
  mu_check(!update);

  update = targets.Trigger(-1000.000);
  mu_check(update);
  mu_check(about_equal_value(targets.value(), -1000));

  targets = ObserverLogTargets(10, 1e-8);
  update = targets.Trigger(1e-9);
  mu_check(update);
  update = targets.Trigger(0);
  mu_check(update);
  update = targets.Trigger(-1.2e-8);
  mu_check(update);
  update = targets.Trigger(-1.2e-7);
  mu_check(update);

  targets = ObserverLogTargets(10, 1e-8);
  update = targets.Trigger(-1.2e-7);
  mu_check(update);
}

/**
 * Tests the function coco_observer_lin_targets_trigger.
 */
MU_TEST(test_coco_observer_lin_targets_trigger) {

  ObserverLinTargets targets(1e-3);
  int update;

  update = targets.Trigger(12);
  mu_check(update);
  mu_check(about_equal_value(targets.value(), 12));

  update = targets.Trigger(10);
  mu_check(update);
  mu_check(about_equal_value(targets.value(), 10));

  update = targets.Trigger(2);
  mu_check(update);
  mu_check(about_equal_value(targets.value(), 2));

  update = targets.Trigger(1.2);
  mu_check(update);
  mu_check(about_equal_value(targets.value(), 1.2));

  update = targets.Trigger(0.12);
  mu_check(update);
  mu_check(about_equal_value(targets.value(), 0.12));

  update = targets.Trigger(10);
  mu_check(!update);

  update = targets.Trigger(2);
  mu_check(!update);

  update = targets.Trigger(0.000013);
  mu_check(update);
  mu_check(about_equal_value(targets.value(), 0.001));

  update = targets.Trigger(0.000012);
  mu_check(!update);

  update = targets.Trigger(1e-8);
  mu_check(!update);

  update = targets.Trigger(0);
  mu_check(update);
  mu_check(about_equal_value(targets.value(), 0));

  update = targets.Trigger(-1.2e-8);
  mu_check(!update);

  update = targets.Trigger(-0.001);
  mu_check(update);
  mu_check(about_equal_value(targets.value(), -0.001));

  update = targets.Trigger(-0.0015);
  mu_check(!update);

}

/**
 * Tests the function coco_observer_evaluations_trigger.
 */
MU_TEST(test_coco_observer_evaluations_trigger) {

  size_t evals[53] = { 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 14, 15, 17, 19, 20, 22, 25, 28, 31, 35, 39, 40, 44,
      50, 56, 63, 70, 79, 89, 100, 112, 125, 141, 158, 177, 199, 200, 223, 251, 281, 316, 354, 398, 400, 446,
      501, 562, 630, 707, 794, 891, 1000};

  size_t i, j;
  int update, found;

  ObserverEvaluations evaluations = ObserverEvaluations(20, "1,2,5", 2); /* do not assume that 20 is default */
  for (i = 1; i <= 1000; i++) {
    update = evaluations.Trigger(i);
    found = 0;
    for (j = 0; j < 53; j++) {
      if (i == evals[j]) {
        found = 1;
        break;
      }
    }
    if (update != found) {
      coco_warning("test_coco_observer_evaluations_trigger(): Assert fails for evaluation number = %lu",
      		(unsigned long) i);
    }
    mu_check(update == found);
  }

}

} // namespace coco

int main(void) {
  MU_RUN_TEST(coco::test_coco_observer_targets_trigger);
  MU_RUN_TEST(coco::test_coco_observer_log_targets_trigger);
  MU_RUN_TEST(coco::test_coco_observer_lin_targets_trigger);
  MU_RUN_TEST(coco::test_coco_observer_evaluations_trigger);

  MU_REPORT();

  return MU_EXIT_CODE;
}
