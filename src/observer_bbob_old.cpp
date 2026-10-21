/**
 * @file observer_bbob_old.cpp
 * @brief Old implementation of the bbob observer.
 */

#include "coco.h"
#include "coco_utilities.cpp"

/**
 * @brief Which info file the bbob_old loggers of one observer write to and what it contains.
 *
 * Shared by the observer and its loggers, which may be freed after the observer.
 */
struct logger_bbob_old_info_state_t {
  /** @brief The number of dimensions an info file can hold. TODO should use dimensions from the suite */
  static constexpr size_t number_of_dimensions = 6;

  size_t current_dimension = 0; /**< @brief Dimension of the last logged problem */
  size_t current_function = 0;  /**< @brief Function of the last logged problem */
  size_t first_instance = 0;    /**< @brief The first instance in the current info file (part of its name) */
  std::string first_instance_string;
  /** @brief The dimensions already in the current info file. If the one about to be logged is new for the
   * current function, a new info file is created with the number of the first instance as suffix. */
  size_t dimensions_in_info_file[number_of_dimensions] = {};
  /**
   * @brief Whether a logger is open. The info file handling fails if several problems are open, so this
   * leads to an error.
   */
  bool logger_is_open = false;
};

/**
 * @brief The bbob_old observer (has no options of its own).
 */
struct observer_bbob_old_t final : public coco_observer_s {
  explicit observer_bbob_old_t(const coco_observer_settings_s &observer_settings)
      : coco_observer_s(observer_settings), info_state(std::make_shared<logger_bbob_old_info_state_t>()) {}

  coco_problem_ptr_t add_logger(coco_problem_ptr_t problem) override;

  std::shared_ptr<logger_bbob_old_info_state_t> info_state; /**< @brief Shared with the loggers */
};
