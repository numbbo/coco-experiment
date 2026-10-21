/**
 * @file observer_bbob.cpp
 * @brief Implementation of the bbob observer.
 */

#include "coco.h"
#include "coco_utilities.cpp"

namespace coco {

class LoggerBbob;

/**
 * @brief The bbob observer.
 *
 * There is a cyclic reference between the bbob_logger and the bbob_observer (through the observer's
 * observed_problem and the logger's data, which points to the observer). This is needed to be able
 * to free both objects without problems.
 *
 * Possible options:
 *
 * - "prefix: STRING" defines the prefix of the name of the info files. The default value is "bbobexp".
 */
class ObserverBbob final : public Observer {
public:
  ObserverBbob(const ObserverSettings &observer_settings, std::string_view options, OptionKeys &option_keys)
      : Observer(observer_settings), prefix_("bbobexp"), last_function_(0), last_dimension_(0), last_dat_file_("init") {

    /* Sets the valid keys for bbob observer options
     * IMPORTANT: This list should be up-to-date with the code and the documentation */
    option_keys = {"prefix"};
    OptionsReadString(options, "prefix", prefix_);
  }

  ProblemPtr AddLogger(ProblemPtr problem) override;

private:
  friend class LoggerBbob; /* Keeps the state of the last logged problem in the observer */

  std::string prefix_;              /**< @brief Prefix in the name of the info and data files */
  /* Store the information on the last function, dimension and dat file used, to be able to tell when the
   * logger needs to write in a new dat file.  */
  size_t last_function_;      /**< @brief The function that was logged last */
  size_t last_dimension_;     /**< @brief The dimension that was logged last */
  std::string last_dat_file_; /**< @brief The name of the .dat file that was last used for logging */
};

} // namespace coco
