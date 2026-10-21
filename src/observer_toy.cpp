/**
 * @file observer_toy.cpp
 * @brief Implementation of the toy observer.
 */

#include "coco.h"
#include "coco_utilities.cpp"

static coco_problem_ptr_t logger_toy(coco_observer_t *observer, coco_problem_ptr_t problem);

/**
 * @brief The toy observer.
 *
 * Possible options:
 * - file_name: string (name of the output file; default value is "first_hitting_times.dat")
 */
struct observer_toy_t final : public coco_observer_s {

  observer_toy_t(const coco_observer_settings_s &observer_settings, std::string_view options,
                 coco_option_keys_t &option_keys)
      : coco_observer_s(observer_settings), log_file(nullptr) {

    std::string file_name = "first_hitting_times.dat";

    /* Sets the valid keys for toy observer options
     * IMPORTANT: This list should be up-to-date with the code and the documentation */
    option_keys = {"file_name"};

    /* Read file_name from the options and use it to initialize the observer */
    coco_options_read_string(options, "file_name", file_name);

    /* Open log_file */
    coco_create_directory(result_folder);
    const std::string path = coco_join_path(result_folder, file_name);

    log_file = fopen(path.c_str(), "a");
    if (log_file == nullptr) {
      coco_error("observer_toy(): failed to open file %s.", path.c_str());
      return; /* Never reached */
    }
  }

  ~observer_toy_t() override {
    if (log_file != nullptr)
      fclose(log_file);
  }

  coco_problem_ptr_t add_logger(coco_problem_ptr_t problem) override { return logger_toy(this, std::move(problem)); }

  FILE *log_file; /**< @brief File used for logging. */
};

/**
 * @brief Returns the toy observer (the given observer must be one).
 */
static observer_toy_t *observer_toy_get(coco_observer_t *observer) {
  assert(dynamic_cast<observer_toy_t *>(observer) != nullptr);
  return static_cast<observer_toy_t *>(observer);
}
