/**
 * @file observer_toy.cpp
 * @brief Implementation of the toy observer.
 */

#include "coco.h"
#include "coco_utilities.cpp"

/**
 * @brief The toy observer.
 *
 * Possible options:
 * - file_name: string (name of the output file; default value is "first_hitting_times.dat")
 */
struct observer_toy_t final : public coco_observer_s {

  observer_toy_t(const coco_observer_settings_s &observer_settings, std::string_view options,
                 coco_option_keys_t &option_keys)
      : coco_observer_s(observer_settings) {

    std::string file_name = "first_hitting_times.dat";

    /* Sets the valid keys for toy observer options
     * IMPORTANT: This list should be up-to-date with the code and the documentation */
    option_keys = {"file_name"};

    /* Read file_name from the options and use it to initialize the observer */
    coco_options_read_string(options, "file_name", file_name);

    /* Open log_file */
    coco_create_directory(result_folder);
    const std::string path = coco_join_path(result_folder, file_name);

    log_file.open(path, "a");
    if (log_file == nullptr) {
      coco_error("observer_toy_t(): failed to open file %s.", path.c_str());
      return; /* Never reached */
    }
  }

  coco_problem_ptr_t add_logger(coco_problem_ptr_t problem) override;

  coco_file_t log_file; /**< @brief File used for logging. */
};

