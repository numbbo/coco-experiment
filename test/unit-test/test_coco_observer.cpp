#include "minunit.h"

#include "coco.cpp"
#include "about_equal.h"

/**
 * Tests the function coco_observer_targets_trigger.
 */
MU_TEST(test_coco_observer_targets_trigger) {

  coco_observer_targets_t targets = coco_observer_targets_t(0, 1e-3, 1, 1e-5);
  int update;

  update = targets.trigger(99.99999);
  mu_check(update);
  mu_check(about_equal_value(targets.last_target(), 100));

  update = targets.trigger(0.99999);
  mu_check(update);
  mu_check(about_equal_value(targets.last_target(), 1));

  update = targets.trigger(0.899999);
  mu_check(update);
  mu_check(about_equal_value(targets.last_target(), 0.9));

  update = targets.trigger(0.0199999);
  mu_check(update);
  mu_check(about_equal_value(targets.last_target(), 0.02));

  update = targets.trigger(0.00000099999);
  mu_check(update);
  mu_check(about_equal_value(targets.last_target(), 1e-3));

  /* An improvement on the lin targets should not trigger an update if it is not an improvement overall */
  update = targets.trigger(0.00999);
  mu_check(!update);

  update = targets.trigger(0);
  mu_check(update);
  mu_check(about_equal_value(targets.last_target(), 0));

  update = targets.trigger(-0.00100001);
  mu_check(update);
  mu_check(about_equal_value(targets.last_target(), -1.0e-3));

  update = targets.trigger(-0.99999);
  mu_check(update);
  mu_check(about_equal_value(targets.last_target(), -0.999));

  update = targets.trigger(-0.100099999);
  mu_check(!update);

}

/**
 * Tests the function coco_observer_log_targets_trigger.
 */
MU_TEST(test_coco_observer_log_targets_trigger) {

  coco_observer_log_targets_t targets = coco_observer_log_targets_t(10, 1e-8);
  int update;

  update = targets.trigger(99.99999);
  mu_check(update);
  mu_check(about_equal_value(targets.value, 100));

  update = targets.trigger(10);
  mu_check(update);
  mu_check(about_equal_value(targets.value, 10));

  update = targets.trigger(1.00000);
  mu_check(update);
  mu_check(about_equal_value(targets.value, 1));

  update = targets.trigger(10);
  mu_check(!update);

  update = targets.trigger(2);
  mu_check(!update);

  update = targets.trigger(0.00001);
  mu_check(update);
  mu_check(about_equal_value(targets.value, 0.00001));

  update = targets.trigger(12);
  mu_check(!update);

  update = targets.trigger(1e-8);
  mu_check(update);
  mu_check(about_equal_value(targets.value,  1e-8));

  update = targets.trigger(1e-9);
  mu_check(!update);

  update = targets.trigger(-1e-8);
  mu_check(update);
  mu_check(about_equal_value(targets.value,  -1e-8));

  update = targets.trigger(-1e-7);
  mu_check(update);
  mu_check(about_equal_value(targets.value,  -1e-7));

  update = targets.trigger(2);
  mu_check(!update);

  update = targets.trigger(-1000.000);
  mu_check(update);
  mu_check(about_equal_value(targets.value,  -1000));


  targets = coco_observer_log_targets_t(10, 1e-8);
  update = targets.trigger(1e-9);
  mu_check(update);
  update = targets.trigger(0);
  mu_check(update);
  update = targets.trigger(-1.2e-8);
  mu_check(update);
  update = targets.trigger(-1.2e-7);
  mu_check(update);

  targets = coco_observer_log_targets_t(10, 1e-8);
  update = targets.trigger(-1.2e-7);
  mu_check(update);
}


/**
 * Tests the function coco_observer_lin_targets_trigger.
 */
MU_TEST(test_coco_observer_lin_targets_trigger) {

  coco_observer_lin_targets_t targets;
  targets.precision = 1e-3;
  int update;

  update = targets.trigger(12);
  mu_check(update);
  mu_check(about_equal_value(targets.value, 12));

  update = targets.trigger(10);
  mu_check(update);
  mu_check(about_equal_value(targets.value, 10));

  update = targets.trigger(2);
  mu_check(update);
  mu_check(about_equal_value(targets.value, 2));

  update = targets.trigger(1.2);
  mu_check(update);
  mu_check(about_equal_value(targets.value, 1.2));

  update = targets.trigger(0.12);
  mu_check(update);
  mu_check(about_equal_value(targets.value, 0.12));

  update = targets.trigger(10);
  mu_check(!update);

  update = targets.trigger(2);
  mu_check(!update);

  update = targets.trigger(0.000013);
  mu_check(update);
  mu_check(about_equal_value(targets.value, 0.001));

  update = targets.trigger(0.000012);
  mu_check(!update);

  update = targets.trigger(1e-8);
  mu_check(!update);

  update = targets.trigger(0);
  mu_check(update);
  mu_check(about_equal_value(targets.value, 0));

  update = targets.trigger(-1.2e-8);
  mu_check(!update);

  update = targets.trigger(-0.001);
  mu_check(update);
  mu_check(about_equal_value(targets.value, -0.001));

  update = targets.trigger(-0.0015);
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

  coco_observer_evaluations_t evaluations = coco_observer_evaluations_t(20, "1,2,5", 2); /* do not assume that 20 is default */
  for (i = 1; i <= 1000; i++) {
    update = evaluations.trigger(i);
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

int main(void) {
  MU_RUN_TEST(test_coco_observer_targets_trigger);
  MU_RUN_TEST(test_coco_observer_log_targets_trigger);
  MU_RUN_TEST(test_coco_observer_lin_targets_trigger);
  MU_RUN_TEST(test_coco_observer_evaluations_trigger);
	
	MU_REPORT();

  return MU_EXIT_CODE;
}
