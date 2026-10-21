/**
 * @file observer_toy.cpp
 * @brief Implementation of the toy observer.
 */

#include "coco.h"
#include "coco_utilities.cpp"

namespace coco {

/**
 * @brief The toy observer.
 *
 * Possible options:
 * - file_name: string (name of the output file; default value is "first_hitting_times.dat")
 */
class ObserverToy final : public Observer {
public:
  ObserverToy(const ObserverSettings &observer_settings, std::string_view options, OptionKeys &option_keys)
      : Observer(observer_settings) {

    std::string file_name = "first_hitting_times.dat";

    /* Sets the valid keys for toy observer options
     * IMPORTANT: This list should be up-to-date with the code and the documentation */
    option_keys = {"file_name"};

    /* Read file_name from the options and use it to initialize the observer */
    OptionsReadString(options, "file_name", file_name);

    /* Open log_file */
    CreateDirectory(result_folder);
    const std::string path = JoinPath(result_folder, file_name);

    log_file_.Open(path, "a");
    if (log_file_ == nullptr) {
      coco_error("ObserverToy(): failed to open file %s.", path.c_str());
      return; /* Never reached */
    }
  }

  ProblemPtr AddLogger(ProblemPtr problem) override;

private:
  friend class LoggerToy; /* Writes to the log file */

  File log_file_; /**< @brief File used for logging. */
};

} // namespace coco
