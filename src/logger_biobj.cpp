/**
 * @file logger_biobj.cpp
 * @brief Implementation of the bbob-biobj logger.
 *
 * Logs the performance of an optimizer on bi-objective problems with or without constraints and with or
 * without knowing their true Pareto front (and set). Uses the hypervolume indicator, but supports
 * inclusion of other indicators. Archives nondominated solutions.
 *
 * In constrained problems only the feasible solutions are logged (except for the first one).
 *
 * Produces these files:
 * - The .info files contain high-level information on the performed experiment.
 * - The .dat files contain function evaluations, indicator values and target hits for every performance
 * target hit as well as for the last evaluation. Logarithmic targets are used only if an estimation of the
 * optimal indicator value is known. If it is not, linear targets are used instead (this is decided on the
 * suite level and passed to the logger during initialization.
 * - The .tdat files contain function evaluation and indicator values for every predefined evaluation
 * number as well as for the last evaluation.
 * - The .rdat files contain function evaluation and indicator values for each restart of the algorithm.
 * - The .mdat file contains function evaluation and indicator values for each recommended solution.
 * - The .adat files are archive files that contain function evaluations, 2 objectives and dim variables
 * for every nondominated solution. Whether these files are created, at what point in time the logger writes
 * nondominated solutions to the archive and whether the decision variables are output or not depends on
 * the values of log_nondom_mode and log_nondom_mode. See the bi-objective observer constructor
 * observer_biobj() for more information.
 *
 * One .info file is created for each problem group (and indicator type) and contains information on all
 * the problems in that problem group (and indicator type).
 * One .dat, .tdat and .rdat file is created for each problem function and dimension (and indicator type)
 * and contains information for all instances of that problem (and indicator type).
 * One .mdat file is created for each problem function and dimension and contains information for all
 * instances of that problem.
 * One .adat file is created for each problem function, dimension and instance.
 *
 * @note Whenever in this file a ROI is mentioned, it means the (normalized) region of interest in the
 * objective space. The non-normalized ROI is a rectangle with the ideal and nadir points as its two
 * opposite vertices, while the normalized ROI is the square [0, 1]^2. If not specifically mentioned, the
 * normalized ROI is assumed.
 */

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "coco.h"
#include "coco_internal.h"

#include "coco_utilities.cpp"
#include "coco_problem.cpp"
#include "coco_string.cpp"
#include "mo_avl_tree.cpp"
#include "observer_biobj.cpp"

#include "mo_utilities.cpp"

/** @brief Number of implemented indicators */
static constexpr size_t LOGGER_BIOBJ_NUMBER_OF_INDICATORS = 1;

/** @brief Names of implemented indicators
 *
 * "hyp" stands for the hypervolume indicator.
 * */
static const std::string_view logger_biobj_indicators[LOGGER_BIOBJ_NUMBER_OF_INDICATORS] = {"hyp"};

/**
 * @brief The indicator type.
 *
 * <B> The hypervolume indicator ("hyp") </B>
 *
 * The hypervolume indicator measures the volume of the portion of the ROI in the objective space that is
 * dominated by the current Pareto front approximation. Instead of logging the hypervolume indicator value,
 * this implementation logs the difference between the best know hypervolume indicator (a value stored in
 * best_value) and the hypervolume indicator of the current Pareto front approximation (current_value). The
 * current_value equals 0 if no solution is located in the ROI. In order to be able to register the
 * performance of an optimizer even before the ROI is reached, an additional value is computed when no
 * solutions are located inside the ROI. This value is stored in additional_penalty and equals the
 * normalized distance to the ROI of the solution closest to the ROI (additional_penalty is set to 0 as
 * soon as a solution reaches the ROI). The final value to be logged (overall_value) is therefore computed
 * in the following way:
 *
 * overall_value = best_value - current_value + additional_penalty
 *
 * If the suite does not provide an estimation for the best hypervolume value, best_value is set to 1.0.
 *
 * @note Other indicators are yet to be implemented.
 */
struct logger_biobj_t;

struct logger_biobj_indicator_t {

  /** @brief Initializes the indicator: opens its files and resets the counters. */
  logger_biobj_indicator_t(const logger_biobj_t &logger, const observer_biobj_t &observer,
                           const coco_problem_t &problem, std::string_view name);

  /** @brief Outputs the final information about this indicator. */
  void finalize(const logger_biobj_t &logger);

  ~logger_biobj_indicator_t() {
    coco_debug("Started logger_biobj_indicator_t::~logger_biobj_indicator_t()");
    if (dat_file != nullptr)
      fclose(dat_file);
    if (tdat_file != nullptr)
      fclose(tdat_file);
    if (rdat_file != nullptr)
      fclose(rdat_file);
    if (info_file != nullptr)
      fclose(info_file);
    coco_debug("Ended   logger_biobj_indicator_t::~logger_biobj_indicator_t()");
  }

  std::string name; /**< @brief Name of the indicator used for identification and the output. */

  FILE *info_file; /**< @brief File for logging summary information on algorithm performance. */
  FILE *dat_file;  /**< @brief File for logging indicator values at predefined values */
  FILE *tdat_file; /**< @brief File for logging indicator values at predefined evaluations. */
  FILE *rdat_file; /**< @brief File for logging restart information */

  bool target_hit; /**< @brief Whether the performance target was hit in the latest evaluation. */
  coco_observer_targets_t targets;
  /**< @brief Triggers based on performance target values. */
  bool evaluation_logged; /**< @brief Whether the whether the latest evaluation was logged. */
  coco_observer_evaluations_t evaluations;
  /**< @brief Triggers based on the number of evaluations. */

  double best_value;         /**< @brief The best known indicator value for this problem. */
  double current_value;      /**< @brief The current indicator value. */
  double additional_penalty; /**< @brief Additional penalty for solutions outside the ROI. */
  double overall_value;      /**< @brief The overall value of the indicator tested for target hits. */
  double previous_value;     /**< @brief The previous overall value of the indicator. */

private:
  /** @brief Creates and initializes one of the *dat files (.dat, .tdat or .rdat). */
  void initialize_file(const logger_biobj_t &logger, const observer_biobj_t &observer, const coco_problem_t &problem,
                       FILE **f, const char *file_ending, const int output_targets);
};

/**
 * @brief The type for the node's item in the AVL tree as used by the bi-objective logger.
 *
 * Contains information on the exact objective values (y) and their rounded normalized values (normalized_y).
 * The exact values are used for output, while archive update and indicator computation use the normalized
 * values.
 */
struct logger_biobj_avl_item_t {
  std::vector<double> x;            /**< @brief The decision values of this solution. */
  std::vector<double> y;            /**< @brief The values of objectives of this solution. */
  bool is_feasible;                 /**< @brief Whether the solution is feasible. */
  std::vector<double> normalized_y; /**< @brief The values of normalized objectives of this solution. */
  size_t evaluation_number;         /**< @brief The evaluation number of when the solution was created. */

  double indicator_contribution[LOGGER_BIOBJ_NUMBER_OF_INDICATORS];
  /**< @brief The contribution of this solution to the overall indicator values. */
  int within_ROI; /**< @brief Whether the solution is within the region of interest (ROI). */
};

/**
 * @brief Creates and returns the information on the solution in the form of a node's item in the AVL tree.
 */
static logger_biobj_avl_item_t *logger_biobj_node_create(const coco_problem_t *problem, const double *x,
                                                         const double *y, const double *constraints,
                                                         const size_t evaluation_number, const size_t dim,
                                                         const size_t num_obj, const size_t num_const) {

  double sum_constraints = 0;

  logger_biobj_avl_item_t *item = new logger_biobj_avl_item_t();

  /* Copy the data */
  item->x.assign(x, x + dim);
  item->y.assign(y, y + num_obj);

  /* Compute the normalized y */
  item->normalized_y = mo_normalize(item->y.data(), problem->best_value.data(), problem->nadir_value.data(), num_obj);
  item->within_ROI = mo_is_within_ROI(item->normalized_y.data(), num_obj);

  item->evaluation_number = evaluation_number;
  for (size_t i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++)
    item->indicator_contribution[i] = 0;

  /* Determine whether the solution is feasible */
  item->is_feasible = 1;
  if (num_const > 0) {
    assert(constraints != nullptr);
    for (size_t i = 0; i < num_const; i++) {
      if (constraints[i] > 0)
        sum_constraints += constraints[i];
    }
    item->is_feasible = !(sum_constraints > 0);
  }

  return item;
}

/**
 * @brief Frees the data of the given logger_biobj_avl_item_t.
 */
static void logger_biobj_node_free(logger_biobj_avl_item_t *item, void *userdata) {

  delete item;
  (void)userdata; /* To silence the compiler */
}

/**
 * @brief The biobj logger, which wraps the observed problem.
 *
 * @note Some fields from the observers (coco_observer as well as observer_biobj) are copied here because the
 * observers can be deleted before the logger is finalized and we need these fields for finalization.
 */
struct logger_biobj_t final : public coco_problem_transformed_t {
  observer_biobj_t *observer{}; /**< @brief Pointer to the observer (might be NULL at the end) */

  observer_biobj_log_nondom_e log_nondom_mode{};
  /**< @brief Mode for archiving nondominated solutions. */
  FILE *adat_file{}; /**< @brief File for archiving nondominated solutions (all or final). */
  FILE *mdat_file{}; /**< @brief File for logging recommended solutions */

  bool log_vars{};           /**< @brief Whether to log the decision values. */
  int precision_x{};         /**< @brief Precision for outputting decision values. */
  int precision_f{};         /**< @brief Precision for outputting objective values. */
  bool log_discrete_as_int{}; /**< @brief Whether to output discrete variables in int or double format. */
  bool algorithm_restarted{}; /**< @brief Whether the algorithm has restarted (output information to .rdat file). */

  size_t num_func_evaluations{}; /**< @brief The number of function evaluations performed so far. */
  size_t num_cons_evaluations{}; /**< @brief The number of evaluations of constraints performed so far. */

  size_t previous_evaluations{}; /**< @brief The number of evaluations from the previous call to the logger. */

  avl_tree_ptr_t archive_tree; /**< @brief The tree keeping currently non-dominated solutions. */
  avl_tree_ptr_t buffer_tree;  /**< @brief The tree with pointers to nondominated solutions that haven't
                                    been logged yet. */

  /* TODO: Implement other indicators */
  bool compute_indicators{}; /**< @brief Whether to compute the indicators. */
  std::unique_ptr<logger_biobj_indicator_t> indicators[LOGGER_BIOBJ_NUMBER_OF_INDICATORS];
  /**< @brief The implemented indicators. */

  logger_biobj_t(observer_biobj_t &observer, coco_problem_ptr_t inner);

  ~logger_biobj_t() override;

  void evaluate_function(const double *x, double *y) override;
  void recommend_solution(const double *x) override;
  int feed_solution(const size_t evaluation, const double *y);

  /** @brief Saves the information that the algorithm has restarted. */
  void signal_restart();

  /** @brief Returns the biobj logger behind the given problem (which must be one). */
  static logger_biobj_t *from(const coco_problem_t *problem);

private:
  /** @brief Adds the solution to the archive if it is not dominated, returns whether the archive changed. */
  int update_archive(logger_biobj_avl_item_t *node_item);

  /** @brief Outputs the indicator information according to the observer options. */
  void output(const int update_performed, const logger_biobj_avl_item_t *node_item);

  /** @brief Outputs the final nondominated solutions to the archive file. */
  void write_final_archive();

  /** @brief The decision values are not transformed. */
  bool inverse(const double *inner_x, double *x) const override { return identity_inverse(inner_x, x); }
};

/**
 * @brief Defines the ordering of AVL tree nodes based on the value of the last objective.
 *
 * @note This ordering is used by the archive_tree.
 */
static int avl_tree_compare_by_last_objective(const logger_biobj_avl_item_t *item1,
                                              const logger_biobj_avl_item_t *item2, void *userdata) {
  (void)userdata; /* To silence the compiler */
  if (coco_double_almost_equal(item1->normalized_y[1], item2->normalized_y[1], mo_precision))
    return 0;
  else if (item1->normalized_y[1] < item2->normalized_y[1])
    return -1;
  else
    return 1;
}

/**
 * @brief Defines the ordering of AVL tree nodes based on the evaluation number (the time when the nodes were
 * created).
 *
 * @note This ordering is used by the buffer_tree.
 */
static int avl_tree_compare_by_eval_number(const logger_biobj_avl_item_t *item1, const logger_biobj_avl_item_t *item2,
                                           void *userdata) {
  (void)userdata; /* To silence the compiler */
  if (item1->evaluation_number < item2->evaluation_number)
    return -1;
  else if (item1->evaluation_number > item2->evaluation_number)
    return 1;
  else
    return 0;
}

/**
 * @brief Outputs the AVL tree to the given file. Returns the number of nodes in the tree.
 */
static size_t logger_biobj_tree_output(FILE *file, const avl_tree_t *tree, const size_t dim, const size_t num_int_vars,
                                       const size_t num_obj, const int log_vars, const int precision_x,
                                       const int precision_f, const int log_discrete_as_int) {

  avl_node_t *solution;
  size_t i;
  size_t number_of_nodes = 0;

  if (tree->tail) {
    /* There is at least a solution in the tree to output */
    solution = tree->head;
    while (solution != nullptr) {
      fprintf(file, "%lu\t", (unsigned long)((logger_biobj_avl_item_t *)solution->item)->evaluation_number);
      for (size_t j = 0; j < num_obj; j++)
        fprintf(file, "%.*e\t", precision_f, ((logger_biobj_avl_item_t *)solution->item)->y[j]);
      if (log_vars) {
        for (i = 0; i < dim; i++)
          if ((i < num_int_vars) && (log_discrete_as_int))
            fprintf(file, "%d\t", coco_double_to_int(((logger_biobj_avl_item_t *)solution->item)->x[i]));
          else
            fprintf(file, "%.*e\t", precision_x, ((logger_biobj_avl_item_t *)solution->item)->x[i]);
      }
      fprintf(file, "\n");
      solution = solution->next;
      number_of_nodes++;
    }
  }

  return number_of_nodes;
}

/**
 * @brief Updates the archive and buffer trees with the given node.
 *
 * Checks for domination and updates the archive tree and the values of the indicators if the given node is
 * not weakly dominated by existing nodes in the archive tree. This is where the main computation of
 * indicator values takes place.
 *
 * @return 1 if the update was performed and 0 otherwise.
 */
int logger_biobj_t::update_archive(logger_biobj_avl_item_t *node_item) {

  avl_node_t *node, *next_node, *new_node;
  int trigger_update = 0;
  int dominance;
  int previous_unavailable = 0;

  /* If the node contains an infeasible solution, exit immediately (do not update the tree) */
  if (node_item->is_feasible == 0)
    return 0;

  /* Find the first point that is not worse than the new point (NULL if such point does not exist) */
  node = avl_item_search_right(archive_tree.get(), node_item, nullptr);

  if (node == nullptr) {
    /* The new point is an extreme point */
    trigger_update = 1;
    next_node = archive_tree->head;
  } else {
    dominance = mo_get_dominance(node_item->normalized_y.data(),
                                 ((logger_biobj_avl_item_t *)node->item)->normalized_y.data(), number_of_objectives);
    if (dominance > -1) {
      trigger_update = 1;
      next_node = node->next;
      if (dominance == 1) {
        /* The new point dominates the next point, remove the next point */
        if (compute_indicators) {
          for (size_t i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
            indicators[i]->current_value -= ((logger_biobj_avl_item_t *)node->item)->indicator_contribution[i];
          }
        }
        avl_item_delete(buffer_tree.get(), node->item);
        avl_node_delete(archive_tree.get(), node);
      }
    } else {
      /* The new point is dominated or equal to an existing one, nothing more to do */
      trigger_update = 0;
    }
  }

  if (!trigger_update) {
    logger_biobj_node_free(node_item, nullptr);
  } else {
    /* Perform tree update */
    while (next_node != nullptr) {
      /* Check the dominance relation between the new node and the next node. There are only two possibilities:
       * dominance = 0: the new node and the next node are nondominated
       * dominance = 1: the new node dominates the next node */
      node = next_node;
      dominance = mo_get_dominance(node_item->normalized_y.data(),
                                   ((logger_biobj_avl_item_t *)node->item)->normalized_y.data(), number_of_objectives);
      if (dominance == 1) {
        /* The new point dominates the next point, remove the next point */
        if (compute_indicators) {
          for (size_t i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
            indicators[i]->current_value -= ((logger_biobj_avl_item_t *)node->item)->indicator_contribution[i];
          }
        }
        next_node = node->next;
        avl_item_delete(buffer_tree.get(), node->item);
        avl_node_delete(archive_tree.get(), node);
      } else {
        break;
      }
    }

    new_node = avl_item_insert(archive_tree.get(), node_item);
    assert(new_node != nullptr);
    avl_item_insert(buffer_tree.get(), node_item);

    if (compute_indicators) {
      if (node_item->within_ROI) {
        /* Compute indicator value for new node and update the indicator value of the affected nodes */
        logger_biobj_avl_item_t *next_item, *previous_item;

        if (new_node->next != nullptr) {
          next_item = (logger_biobj_avl_item_t *)new_node->next->item;
          if (next_item->within_ROI) {
            for (size_t i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
              indicators[i]->current_value -= next_item->indicator_contribution[i];
              if (indicators[i]->name == "hyp") {
                next_item->indicator_contribution[i] =
                    (node_item->normalized_y[0] - next_item->normalized_y[0]) * (1 - next_item->normalized_y[1]);
                assert(next_item->indicator_contribution[i] >= 0);
              } else {
                coco_error(
                    "logger_biobj_t::update_archive(): Indicator computation not implemented yet for indicator %s",
                    indicators[i]->name.c_str());
              }
              indicators[i]->current_value += next_item->indicator_contribution[i];
            }
          }
        }

        previous_unavailable = 0;
        if (new_node->prev != nullptr) {
          previous_item = (logger_biobj_avl_item_t *)new_node->prev->item;
          if (previous_item->within_ROI) {
            for (size_t i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
              if (indicators[i]->name == "hyp") {
                node_item->indicator_contribution[i] =
                    (previous_item->normalized_y[0] - node_item->normalized_y[0]) * (1 - node_item->normalized_y[1]);
                assert(node_item->indicator_contribution[i] >= 0);
              } else {
                coco_error(
                    "logger_biobj_t::update_archive(): Indicator computation not implemented yet for indicator %s",
                    indicators[i]->name.c_str());
              }
            }
          } else {
            previous_unavailable = 1;
          }
        } else {
          previous_unavailable = 1;
        }

        if (previous_unavailable) {
          /* Previous item does not exist or is out of ROI, use reference point instead */
          for (size_t i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
            if (indicators[i]->name == "hyp") {
              node_item->indicator_contribution[i] =
                  (1 - node_item->normalized_y[0]) * (1 - node_item->normalized_y[1]);
              assert(node_item->indicator_contribution[i] >= 0);
            } else {
              coco_error("logger_biobj_t::update_archive(): Indicator computation not implemented yet for indicator %s",
                         indicators[i]->name.c_str());
            }
          }
        }

        for (size_t i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
          if (indicators[i]->name == "hyp") {
            assert(node_item->indicator_contribution[i] >= 0);
            indicators[i]->current_value += node_item->indicator_contribution[i];
          }
        }
      }
    }
  }

  return trigger_update;
}

/**
 * @brief Creates and initializes one of the *dat files (.dat, .tdat or .rdat).
 */
void logger_biobj_indicator_t::initialize_file(const logger_biobj_t &logger, const observer_biobj_t &observer,
                                               const coco_problem_t &problem, FILE **f, const char *file_ending,
                                               const int output_targets) {
  static const char *header_w_targets = "%%\n"
                                        "%% index = %lu, name = %s\n"
                                        "%% instance = %lu, reference value = %.*e\n"
                                        "%% function evaluation | indicator value | target hit\n";
  static const char *header_wo_targets = "%%\n"
                                         "%% index = %lu, name = %s\n"
                                         "%% instance = %lu, reference value = %.*e\n"
                                         "%% function evaluation | indicator value\n";

  /* Create and open the file */
  const std::string folder_name = coco_join_path(observer.result_folder, problem.problem_type);
  coco_create_directory(folder_name);
  const std::string prefix = coco_remove_from_string(problem.problem_id, "_i", "_d");
  const std::string path_name =
      coco_join_path(folder_name, coco_string_format("%s_%s.%s", prefix.c_str(), name.c_str(), file_ending));
  *f = fopen(path_name.c_str(), "a");
  if (*f == nullptr) {
    coco_error("logger_biobj_indicator_t::initialize_file() failed to open file '%s'.", path_name.c_str());
  }

  /* Output header */
  if (output_targets)
    fprintf(*f, header_w_targets, (unsigned long)problem.suite_dep_index, problem.problem_name.c_str(),
            (unsigned long)problem.suite_dep_instance, logger.precision_f, best_value);
  else
    fprintf(*f, header_wo_targets, (unsigned long)problem.suite_dep_index, problem.problem_name.c_str(),
            (unsigned long)problem.suite_dep_instance, logger.precision_f, best_value);
}

/**
 * @brief Initializes the indicator with name indicator_name.
 *
 * Opens files for writing and resets counters.
 */
logger_biobj_indicator_t::logger_biobj_indicator_t(const logger_biobj_t &logger, const observer_biobj_t &observer,
                                                   const coco_problem_t &problem, std::string_view name)
    : name(name), info_file(nullptr), dat_file(nullptr), tdat_file(nullptr), rdat_file(nullptr) {

  const coco_suite_t *suite;
  int info_file_exists = 0;

  coco_debug("Started logger_biobj_indicator_t()");

  assert(problem.suite != nullptr);
  suite = problem.suite;
  target_hit = 0;
  evaluation_logged = 0;
  current_value = 0;
  additional_penalty = DBL_MAX;
  overall_value = 0;
  previous_value = 0;

  /* The best known value of the indicator */
  if (this->name == "hyp")
    suite_biobj_get_best_hyp_value(suite->known_optima, problem.problem_id, &best_value);
  else
    coco_error("logger_biobj_indicator_t(): indicator %s not supported", this->name.c_str());

  targets = coco_observer_targets_t(suite->known_optima, observer.lin_target_precision, observer.number_target_triggers,
                                    observer.log_target_precision);
  evaluations = coco_observer_evaluations_t(observer.number_evaluation_triggers, observer.base_evaluation_triggers,
                                            problem.number_of_variables);

  /* Prepare the info file */
  coco_create_directory(observer.result_folder);
  const std::string path_name = coco_join_path(
      observer.result_folder, coco_string_format("%s_%s.info", problem.problem_type.c_str(), this->name.c_str()));
  info_file_exists = coco_file_exists(path_name.c_str());
  info_file = fopen(path_name.c_str(), "a");
  if (info_file == nullptr) {
    coco_error("logger_biobj_indicator_t() failed to open file '%s'.", path_name.c_str());
    return; /* Never reached */
  }

  /* Output header information to the info file */
  const std::string prefix = coco_remove_from_string(problem.problem_id, "_i", "_d");
  if (!info_file_exists) {
    /* Output algorithm name */
    assert(problem.suite != nullptr);
    fprintf(info_file,
            "suite = '%s', algorithm = '%s', indicator = '%s', folder = '%s', coco_version = '%s', "
            "settings = '%s' \n%% %s",
            problem.suite->suite_name.c_str(), observer.algorithm_name.c_str(), this->name.c_str(),
            problem.problem_type.c_str(), coco_version, observer.settings.c_str(), observer.algorithm_info.c_str());
    if (logger.log_nondom_mode == LOG_NONDOM_READ)
      fprintf(info_file, " (reconstructed)");
  }
  if ((observer.previous_function != (long)problem.suite_dep_function) ||
      (observer.previous_dimension != (long)problem.number_of_variables)) {
    fprintf(info_file, "\nfunction = %2lu, ", (unsigned long)problem.suite_dep_function);
    fprintf(info_file, "dim = %2lu, ", (unsigned long)problem.number_of_variables);
    fprintf(info_file, "%s_%s.dat", prefix.c_str(), this->name.c_str());
  }

  /* Initialize the .dat, .tdat and .rdat files */
  initialize_file(logger, observer, problem, &dat_file, "dat", 1);
  initialize_file(logger, observer, problem, &tdat_file, "tdat", 0);
  initialize_file(logger, observer, problem, &rdat_file, "rdat", 0);

  coco_debug("Ended   logger_biobj_indicator_t()");
}

/**
 * @brief Outputs the final information about this indicator.
 */
void logger_biobj_indicator_t::finalize(const logger_biobj_t &logger) {

  /* Log the last eval_number in the dat file if wasn't already logged */
  if (!target_hit) {
    fprintf(dat_file, "%lu\t%.*e\t%.*e\n", (unsigned long)logger.num_func_evaluations, logger.precision_f,
            overall_value, logger.precision_f, targets.last_target());
  }

  /* Log the last eval_number in the tdat file if wasn't already logged */
  if (!evaluation_logged) {
    fprintf(tdat_file, "%lu\t%.*e\n", (unsigned long)logger.num_func_evaluations, logger.precision_f, overall_value);
  }

  /* Log the information in the info file */
  fprintf(info_file, ", %lu:%lu|%.1e", (unsigned long)logger.suite_dep_instance,
          (unsigned long)logger.num_func_evaluations, overall_value);
  fflush(info_file);
}

/*
 * @brief Outputs the information according to the observer options.
 *
 * Outputs to the:
 * - dat file, if the archive was updated and a new target was reached for an indicator;
 * - tdat file, if the number of evaluations matches one of the predefined numbers.
 *
 * Note that a target is reached when
 * best_value - current_value + additional_penalty <= relative_target_value
 *
 * The relative_target_value is a target for indicator difference, not the actual indicator value!
 */
void logger_biobj_t::output(const int update_performed, const logger_biobj_avl_item_t *node_item) {

  logger_biobj_indicator_t *indicator;

  coco_debug("Started logger_biobj_t::output()");

  if (compute_indicators) {
    for (size_t i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {

      indicator = indicators[i].get();
      indicator->target_hit = 0;
      indicator->previous_value = indicator->overall_value;

      /* If the update was performed, update the overall indicator value */
      if (update_performed) {
        /* Compute the overall_value of an indicator */
        if (indicator->name == "hyp") {
          if (coco_double_almost_equal(indicator->current_value, 0, mo_precision)) {
            /* Update the additional penalty for hypervolume (the minimal distance from the nondominated set
             * to the ROI) */
            double new_distance = mo_get_distance_to_ROI(node_item->normalized_y.data(), number_of_objectives);
            indicator->additional_penalty = coco_double_min(indicator->additional_penalty, new_distance);
            assert(indicator->additional_penalty >= 0);
          } else {
            indicator->additional_penalty = 0;
          }
          indicator->overall_value = indicator->best_value - indicator->current_value + indicator->additional_penalty;
        } else {
          coco_error("logger_biobj_t::output(): Indicator computation not implemented yet for indicator %s",
                     indicator->name.c_str());
        }

        /* Check whether a target was hit */
        indicator->target_hit = indicator->targets.trigger(indicator->overall_value);
      } else if ((num_func_evaluations == 1) && (node_item->is_feasible == 0)) {
        /* Special case if the solution is infeasible and this is the first evaluation */
        indicator->overall_value = INFINITY_FOR_LOGGING;
        indicator->target_hit = indicator->targets.trigger(indicator->overall_value);
      }

      /* Log to the dat file if a performance target was hit */
      if (indicator->target_hit) {
        fprintf(indicator->dat_file, "%lu\t%.*e\t%.*e\n", (unsigned long)num_func_evaluations, precision_f,
                indicator->overall_value, precision_f, indicator->targets.last_target());
      }

      if (log_nondom_mode == LOG_NONDOM_READ) {
        /* Log to the tdat file the previous indicator value if any evaluation number between the previous and
         * this one matches one of the predefined evaluation numbers. */
        for (size_t j = previous_evaluations + 1; j < num_func_evaluations; j++) {
          indicator->evaluation_logged = indicator->evaluations.trigger(j);
          if (indicator->evaluation_logged) {
            fprintf(indicator->tdat_file, "%lu\t%.*e\n", (unsigned long)j, precision_f, indicator->previous_value);
          }
        }
      }

      /* Log to the tdat file if the number of evaluations matches one of the predefined numbers */
      indicator->evaluation_logged = indicator->evaluations.trigger(num_func_evaluations);
      if (indicator->evaluation_logged) {
        fprintf(indicator->tdat_file, "%lu\t%.*e\n", (unsigned long)num_func_evaluations, precision_f,
                indicator->overall_value);
      }

      /* Log to the rdat file if the algorithm was restarted */
      if (algorithm_restarted) {
        fprintf(indicator->rdat_file, "%lu\t%.*e\n", (unsigned long)num_func_evaluations, precision_f,
                indicator->overall_value);
      }
    }
  }

  if (algorithm_restarted)
    algorithm_restarted = 0;

  coco_debug("Ended   logger_biobj_t::output()");
}

/**
 * @brief Evaluates the function, increases the number of evaluations and outputs information according to
 * observer options.
 */
void logger_biobj_t::evaluate_function(const double *x, double *y) {
  logger_biobj_avl_item_t *node_item;
  int update_performed;

  /* Evaluate the objectives */
  coco_evaluate_function(inner_problem.get(), x, y);
  num_func_evaluations++;

  /* Evaluate the constraints */
  std::vector<double> constraints(number_of_constraints);
  if (number_of_constraints > 0)
    inner_problem->evaluate_constraint(x, constraints.data(), 0);
  num_cons_evaluations = evaluations_constraints;

  node_item = logger_biobj_node_create(inner_problem.get(), x, y, constraints.data(), num_func_evaluations,
                                       number_of_variables, number_of_objectives, number_of_constraints);

  /* Update the archive with the new solution, if it is not dominated by or equal to existing solutions in
   * the archive */
  update_performed = update_archive(node_item);

  /* If the archive was updated and you need to log all nondominated solutions, output the new solution to
   * nondom_file */
  if (update_performed && (log_nondom_mode == LOG_NONDOM_ALL)) {
    logger_biobj_tree_output(adat_file, buffer_tree.get(), number_of_variables, number_of_integer_variables,
                             number_of_objectives, log_vars, precision_x, precision_f, log_discrete_as_int);
    avl_tree_purge(buffer_tree.get());

    /* Flush output so that impatient users can see progress. */
    fflush(adat_file);
  }

  /* Output according to observer options */
  output(update_performed, node_item);
}

/**
 * @brief Evaluates the function and outputs information according to observer options to the file with
 * recommendations. The evaluation result is not returned and the evaluation counter is not increased.
 */
void logger_biobj_t::recommend_solution(const double *x) {
  size_t i;

  /* Evaluate the objectives */
  std::vector<double> y(number_of_objectives);
  coco_evaluate_function(inner_problem.get(), x, y.data());

  /* Evaluate the constraints */
  std::vector<double> constraints(number_of_constraints);
  if (number_of_constraints > 0)
    inner_problem->evaluate_constraint(x, constraints.data(), 0);
  num_cons_evaluations = evaluations_constraints;

  /* If this is the first call to the recommend function (on any problem), warn about the large
  resulting files (the observer remembers the warning, unless it has already been freed) */
  if ((observer != nullptr) && (observer->large_mdat_file_warning == 0)) {
    coco_warning("Using the logger_biobj_recommend function may result in the creation of large "
                 "mdat files since every recommended solution is logged");
    observer->large_mdat_file_warning = 1;
  }

  /* Log to the mdat file */
  fprintf(mdat_file, "%lu\t", (unsigned long)num_func_evaluations);
  for (size_t j = 0; j < number_of_objectives; j++)
    fprintf(mdat_file, "%.*e\t", precision_f, y[j]);
  for (size_t j = 0; j < number_of_constraints; j++)
    fprintf(mdat_file, "%.*e\t", precision_f, constraints[j]);
  if (log_vars) {
    for (i = 0; i < number_of_variables; i++)
      if ((i < number_of_integer_variables) && (log_discrete_as_int))
        fprintf(mdat_file, "%d\t", coco_double_to_int(x[i]));
      else
        fprintf(mdat_file, "%.*e\t", precision_x, x[i]);
  }
  fprintf(mdat_file, "\n");
}

/**
 * Sets the number of evaluations, adds the objective vector to the archive and outputs information according
 * to observer options (but does not output the archive).
 *
 * @note Vector y must point to a correctly sized allocated memory region and the given evaluation number must
 * be larger than the existing one.
 *
 * Constraints are not supported.
 *
 * @param evaluation The number of evaluations.
 * @param y The objective vector.
 * @return 1 if archive was updated was done and 0 otherwise.
 */
int logger_biobj_t::feed_solution(const size_t evaluation, const double *y) {
  logger_biobj_avl_item_t *node_item;
  int update_performed;

  assert(log_nondom_mode == LOG_NONDOM_READ);

  /* Set the number of evaluations */
  previous_evaluations = num_func_evaluations;
  if (previous_evaluations >= evaluation)
    coco_error("logger_biobj_t::feed_solution(): Evaluation %lu came before evaluation %lu. Note that "
               "the evaluations need to be always increasing.",
               previous_evaluations, evaluation);
  num_func_evaluations = evaluation;

  /* Update the archive with the new solution */
  std::vector<double> x(number_of_variables);
  for (size_t i = 0; i < number_of_variables; i++)
    x[i] = 0;
  node_item = logger_biobj_node_create(inner_problem.get(), x.data(), y, nullptr, num_func_evaluations,
                                       number_of_variables, number_of_objectives, 0);

  /* Update the archive */
  update_performed = update_archive(node_item);

  /* Output according to observer options */
  output(update_performed, node_item);

  return update_performed;
}

/**
 * @brief Outputs the final nondominated solutions to the archive file.
 */
void logger_biobj_t::write_final_archive() {

  avl_tree_t *resorted_tree;
  avl_node_t *solution;

  coco_debug("Started logger_biobj_t::write_final_archive()");

  /* Re-sort archive_tree according to time stamp and then output it */
  resorted_tree = avl_tree_construct((avl_compare_t)avl_tree_compare_by_eval_number, nullptr);

  if (archive_tree->tail) {
    /* There is at least a solution in the tree to output */
    solution = archive_tree->head;
    while (solution != nullptr) {
      avl_item_insert(resorted_tree, solution->item);
      solution = solution->next;
    }
  }

  logger_biobj_tree_output(adat_file, resorted_tree, number_of_variables, number_of_integer_variables,
                           number_of_objectives, log_vars, precision_x, precision_f, log_discrete_as_int);

  avl_tree_destruct(resorted_tree);

  coco_debug("Ended   logger_biobj_t::write_final_archive()");
}

/**
 * @brief Outputs the final information and closes the files.
 */
logger_biobj_t::~logger_biobj_t() {

  coco_debug("Started logger_biobj_t::~logger_biobj_t()");

  if (log_nondom_mode == LOG_NONDOM_FINAL) {
    write_final_archive();
  }

  if (compute_indicators) {
    for (size_t i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
      indicators[i]->finalize(*this);
      indicators[i].reset();
    }
  }

  if (((log_nondom_mode == LOG_NONDOM_ALL) || (log_nondom_mode == LOG_NONDOM_FINAL)) && (adat_file != nullptr)) {
    fprintf(adat_file, "%% evaluations = %lu\n", (unsigned long)num_func_evaluations);
    fclose(adat_file);
    adat_file = nullptr;
  }

  if (mdat_file != nullptr) {
    fclose(mdat_file);
    mdat_file = nullptr;
  }

  /* The observer might not exist at this point */
  if ((observer != nullptr) && (observer->is_active == 1)) {
    /* If the observer still exists (if it does not, observed_problem does not matter any longer) */
    observer->observed_problem = nullptr;
  }

  coco_debug("Ended   logger_biobj_t::~logger_biobj_t()");
}

/**
 * @brief Saves the information that the algorithm has restarted
 */
void logger_biobj_t::signal_restart() {
  if (num_func_evaluations > 0)
    algorithm_restarted = 1;
}

/**
 * @brief Returns the biobj logger behind the given problem (which must be one).
 */
logger_biobj_t *logger_biobj_t::from(const coco_problem_t *problem) {
  assert(dynamic_cast<const logger_biobj_t *>(problem) != nullptr);
  return static_cast<logger_biobj_t *>(const_cast<coco_problem_t *>(problem));
}

/**
 * @brief Feeds the solution to the bi-objective logger of the given problem, see
 * logger_biobj_t::feed_solution().
 */
int coco_logger_biobj_feed_solution(coco_problem_t *problem, const size_t evaluation, const double *y) {
  return logger_biobj_t::from(problem)->feed_solution(evaluation, y);
}

void observer_biobj_t::signal_restart(coco_problem_t *problem) { logger_biobj_t::from(problem)->signal_restart(); }

/**
 * @brief Makes the observed logger (if any) forget this observer.
 */
observer_biobj_t::~observer_biobj_t() {
  coco_debug("Started observer_biobj_t::~observer_biobj_t()");
  /* Make sure that the observed problem's pointer to the observer points to NULL */
  if (observed_problem != nullptr) {
    observed_problem->observer = nullptr;
    observed_problem = nullptr;
  }
  coco_debug("Ended   observer_biobj_t::~observer_biobj_t()");
}

/**
 * @brief Initializes the biobjective logger.
 *
 * Copies all observer field values that are needed after initialization into logger field values
 * because if the observer is deleted before the suite, the observer is not available anymore when
 * the logger is finalized.
 */
logger_biobj_t::logger_biobj_t(observer_biobj_t &observer, coco_problem_ptr_t inner)
    : coco_problem_transformed_t(std::move(inner), observer.observer_name) {
  const char nondom_folder_name[] = "archive";

  coco_debug("Started logger_biobj_t()");

  if (inner_problem->number_of_objectives != 2) {
    coco_error("logger_biobj_t(): The bi-objective logger cannot log a problem with %lu objective(s)",
               (unsigned long)inner_problem->number_of_objectives);
    return; /* Never reached */
  }

  if (observer.observed_problem != nullptr) {
    coco_error("logger_biobj_t(): The observed problem must be closed before a new problem can be observed");
  }

  this->observer = &observer;
  num_func_evaluations = 0;
  num_cons_evaluations = 0;
  previous_evaluations = 0;
  algorithm_restarted = 0;

  /* Copy values from the observer that you might need even if it does not exist any more */
  log_nondom_mode = observer.log_nondom_mode;
  compute_indicators = observer.compute_indicators;
  precision_x = observer.precision_x;
  precision_f = observer.precision_f;
  log_discrete_as_int = observer.log_discrete_as_int;

  if (((observer.log_vars_mode == LOG_VARS_LOW_DIM) && (inner_problem->number_of_variables > 5)) ||
      (observer.log_vars_mode == LOG_VARS_NEVER))
    log_vars = 0;
  else
    log_vars = 1;

  /* Initialize logging of nondominated solutions into the archive file */
  if ((log_nondom_mode == LOG_NONDOM_ALL) || (log_nondom_mode == LOG_NONDOM_FINAL)) {

    /* Create the path to the file */
    const std::string folder_name = coco_join_path(observer.result_folder, nondom_folder_name);
    coco_create_directory(folder_name);

    /* Construct file name */
    const std::string path_name = coco_join_path(
        folder_name,
        inner_problem->problem_id + ((log_nondom_mode == LOG_NONDOM_ALL) ? "_nondom_all.adat" : "_nondom_final.adat"));

    /* Open and initialize the archive file */
    adat_file = fopen(path_name.c_str(), "a");
    if (adat_file == nullptr) {
      coco_error("logger_biobj_t() failed to open file '%s'.", path_name.c_str());
      return; /* Never reached */
    }

    /* Output header information */
    fprintf(adat_file, "%% instance = %lu, name = %s\n", (unsigned long)inner_problem->suite_dep_instance,
            inner_problem->problem_name.c_str());
    if (log_vars) {
      fprintf(adat_file, "%% function evaluation | %lu objectives | %lu variables\n",
              (unsigned long)inner_problem->number_of_objectives, (unsigned long)inner_problem->number_of_variables);
    } else {
      fprintf(adat_file, "%% function evaluation | %lu objectives \n",
              (unsigned long)inner_problem->number_of_objectives);
    }
  }

  /* Initialize the AVL trees */
  archive_tree.reset(
      avl_tree_construct((avl_compare_t)avl_tree_compare_by_last_objective, (avl_free_t)logger_biobj_node_free));
  buffer_tree.reset(avl_tree_construct((avl_compare_t)avl_tree_compare_by_eval_number, nullptr));

  /* Initialize the indicators */
  if (compute_indicators) {
    for (size_t i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++)
      indicators[i] =
          std::make_unique<logger_biobj_indicator_t>(*this, observer, *inner_problem, logger_biobj_indicators[i]);

    observer.previous_function = (long)inner_problem->suite_dep_function;
    observer.previous_dimension = (long)inner_problem->number_of_variables;
  }

  /* Initialize the file with recommendations */

  /* Create the path to the file */
  const std::string folder_name = coco_join_path(observer.result_folder, problem_type);
  coco_create_directory(folder_name);

  /* Construct file name */
  const std::string prefix = coco_remove_from_string(problem_id, "_i", "_d");
  const std::string path_name = coco_join_path(folder_name, prefix + ".mdat");

  /* Open and initialize the recommendation file */
  mdat_file = fopen(path_name.c_str(), "a");
  if (mdat_file == nullptr) {
    coco_error("logger_biobj_t() failed to open file '%s'.", path_name.c_str());
    return; /* Never reached */
  }

  /* Output header information */
  fprintf(mdat_file, "%% instance = %lu, name = %s\n%% function evaluation | %lu objectives",
          (unsigned long)inner_problem->suite_dep_instance, inner_problem->problem_name.c_str(),
          (unsigned long)inner_problem->number_of_objectives);
  if (inner_problem->number_of_constraints > 0)
    fprintf(mdat_file, " | %lu constraints", (unsigned long)inner_problem->number_of_constraints);
  if (log_vars)
    fprintf(mdat_file, " | %lu variables", (unsigned long)inner_problem->number_of_variables);
  fprintf(mdat_file, "\n");

  coco_debug("Ended   logger_biobj_t()");

  observer.observed_problem = this;
}

coco_problem_ptr_t observer_biobj_t::add_logger(coco_problem_ptr_t problem) {
  return std::make_unique<logger_biobj_t>(*this, std::move(problem));
}
