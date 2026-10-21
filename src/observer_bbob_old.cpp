/**
 * @file observer_bbob_old.cpp
 * @brief Old implementation of the bbob observer.
 */

#include "coco.h"
#include "coco_utilities.cpp"

namespace coco {

/**
 * @brief Which info file the bbob_old loggers of one observer write to and what it contains.
 *
 * Shared by the observer and its loggers, which may be freed after the observer.
 */
struct LoggerBbobOldInfoState {
  /** @brief The number of dimensions an info file can hold. TODO should use dimensions from the suite */
  static constexpr size_t NUMBER_OF_DIMENSIONS = 6;

  size_t current_dimension = 0; /**< @brief Dimension of the last logged problem */
  size_t current_function = 0;  /**< @brief Function of the last logged problem */
  size_t first_instance = 0;    /**< @brief The first instance in the current info file (part of its name) */
  std::string first_instance_string;
  /** @brief The dimensions already in the current info file. If the one about to be logged is new for the
   * current function, a new info file is created with the number of the first instance as suffix. */
  size_t dimensions_in_info_file[NUMBER_OF_DIMENSIONS] = {};
};

/**
 * @brief The bbob_old observer (has no options of its own).
 */
class ObserverBbobOld final : public Observer {
public:
  explicit ObserverBbobOld(const ObserverSettings &observer_settings)
      : Observer(observer_settings), info_state_(std::make_shared<LoggerBbobOldInfoState>()) {}

  ProblemPtr AddLogger(ProblemPtr problem) override;

private:
  friend class LoggerBbobOld; /* Shares the info file state with the observer */

  std::shared_ptr<LoggerBbobOldInfoState> info_state_; /**< @brief Shared with the loggers */
};

} // namespace coco
