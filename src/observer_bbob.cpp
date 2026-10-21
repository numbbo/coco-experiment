/**
 * @file observer_bbob.cpp
 * @brief Implementation of the bbob observer.
 */

#include "coco.h"
#include "coco_utilities.cpp"

struct logger_bbob_t;

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
struct observer_bbob_t final : public coco_observer_s {

  observer_bbob_t(const coco_observer_settings_s &observer_settings, std::string_view options,
                  coco_option_keys_t &option_keys)
      : coco_observer_s(observer_settings), observed_problem(nullptr), prefix("bbobexp"), last_function(0),
        last_dimension(0), last_dat_file("init") {

    /* Sets the valid keys for bbob observer options
     * IMPORTANT: This list should be up-to-date with the code and the documentation */
    option_keys = {"prefix"};
    coco_options_read_string(options, "prefix", prefix);
  }

  ~observer_bbob_t() override;

  coco_problem_ptr_t add_logger(coco_problem_ptr_t problem) override;

  void signal_restart(coco_problem_t *problem) override;

  logger_bbob_t *observed_problem;  /**< @brief The logger of the observed problem (nullptr if none is observed) */
  std::string prefix;               /**< @brief Prefix in the name of the info and data files */

  /* Store the information on the last function, dimension and dat file used, to be able to tell when the
   * logger needs to write in a new dat file.  */
  size_t last_function;      /**< @brief The function that was logged last */
  size_t last_dimension;     /**< @brief The dimension that was logged last */
  std::string last_dat_file; /**< @brief The name of the .dat file that was last used for logging */
};

