/**
 * @file observer_bbob_old.cpp
 * @brief Old implementation of the bbob observer.
 */

#include "coco.h"
#include "coco_utilities.cpp"

static coco_problem_t *logger_bbob_old(coco_observer_t *observer, coco_problem_t *problem);

/**
 * @brief The bbob_old observer (has no options of its own).
 */
struct observer_bbob_old_t final : public coco_observer_s {
  explicit observer_bbob_old_t(const coco_observer_settings_s &observer_settings)
      : coco_observer_s(observer_settings) {}

  coco_problem_t *add_logger(coco_problem_t *problem) override { return logger_bbob_old(this, problem); }
};
