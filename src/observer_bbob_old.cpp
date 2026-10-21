/**
 * @file observer_bbob_old.cpp
 * @brief Old implementation of the bbob observer.
 */

#include "coco.h"
#include "coco_utilities.cpp"


/**
 * @brief The bbob_old observer (has no options of its own).
 */
struct observer_bbob_old_t final : public coco_observer_s {
  explicit observer_bbob_old_t(const coco_observer_settings_s &observer_settings)
      : coco_observer_s(observer_settings) {}

  coco_problem_ptr_t add_logger(coco_problem_ptr_t problem) override;
};
