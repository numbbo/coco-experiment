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

namespace coco {

/** @brief Number of implemented indicators */
static constexpr size_t LOGGER_BIOBJ_NUMBER_OF_INDICATORS = 1;

/** @brief Names of implemented indicators
 *
 * "hyp" stands for the hypervolume indicator.
 * */
static const std::string_view LOGGER_BIOBJ_INDICATORS[LOGGER_BIOBJ_NUMBER_OF_INDICATORS] = {"hyp"};

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
class LoggerBiobj;
struct LoggerBiobjAvlItem;

class LoggerBiobjIndicator {
public:
  /** @brief Initializes the indicator: opens its files and resets the counters. */
  LoggerBiobjIndicator(const LoggerBiobj &logger, const ObserverBiobj &observer, const Problem &problem,
                       std::string_view name);

  /** @brief Outputs the final information about this indicator. */
  void Finalize(const LoggerBiobj &logger);

  /**
   * @brief Updates the overall value of the indicator after an evaluation and outputs it according to the
   * observer options (see LoggerBiobj::output()).
   */
  void Output(const LoggerBiobj &logger, const int update_performed, const LoggerBiobjAvlItem *node_item);

  /** @brief Adds the contribution of a solution that entered the archive to the current value. */
  void AddContribution(const double contribution) { current_value_ += contribution; }

  /** @brief Removes the contribution of a solution that left the archive from the current value. */
  void RemoveContribution(const double contribution) { current_value_ -= contribution; }

  const std::string &name() const { return name_; }
  double current_value() const { return current_value_; }
  double overall_value() const { return overall_value_; }

  /** @brief Sets the best known indicator value (the reference value), which the tests use. */
  void set_best_value(const double value) { best_value_ = value; }

private:
  /** @brief Creates and initializes one of the *dat files (.dat, .tdat or .rdat). */
  void InitializeFile(const LoggerBiobj &logger, const ObserverBiobj &observer, const Problem &problem, File &f,
                      const char *file_ending, const int output_targets);

  std::string name_;      /**< @brief Name of the indicator used for identification and the output. */
  File info_file_;        /**< @brief File for logging summary information on algorithm performance. */
  File dat_file_;         /**< @brief File for logging indicator values at predefined values */
  File tdat_file_;        /**< @brief File for logging indicator values at predefined evaluations. */
  File rdat_file_;        /**< @brief File for logging restart information */
  bool target_hit_;       /**< @brief Whether the performance target was hit in the latest evaluation. */
  ObserverTargets targets_;
  /**< @brief Triggers based on performance target values. */
  bool evaluation_logged_; /**< @brief Whether the whether the latest evaluation was logged. */
  ObserverEvaluations evaluations_;
  double best_value_;         /**< @brief The best known indicator value for this problem. */
  double current_value_;      /**< @brief The current indicator value. */
  double additional_penalty_; /**< @brief Additional penalty for solutions outside the ROI. */
  double overall_value_;      /**< @brief The overall value of the indicator tested for target hits. */
  double previous_value_;     /**< @brief The previous overall value of the indicator. */
};

/**
 * @brief The type for the node's item in the AVL tree as used by the bi-objective logger.
 *
 * Contains information on the exact objective values (y) and their rounded normalized values (normalized_y).
 * The exact values are used for output, while archive update and indicator computation use the normalized
 * values.
 */
struct LoggerBiobjAvlItem {
  /**
   * @brief Stores the solution with its objective values normalized with the ideal and nadir point of the
   * problem and whether it is feasible w.r.t. the given constraint values.
   */
  LoggerBiobjAvlItem(const Problem &problem, std::span<const double> x, std::span<const double> y,
                     std::span<const double> constraints, const size_t evaluation_number)
      : x(x.begin(), x.end()), y(y.begin(), y.end()), is_feasible(true),
        normalized_y(MoNormalize(this->y, problem.best_value(), problem.nadir_value())),
        evaluation_number(evaluation_number), indicator_contribution(), within_ROI(MoIsWithinRoi(normalized_y)) {
    /* Determine whether the solution is feasible */
    if (!constraints.empty()) {
      double sum_constraints = 0;
      for (const double constraint : constraints) {
        if (constraint > 0)
          sum_constraints += constraint;
      }
      is_feasible = !(sum_constraints > 0);
    }
  }

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
 * @brief Frees the data of the given LoggerBiobjAvlItem.
 */
static void LoggerBiobjNodeFree(void *item, void *) { delete static_cast<LoggerBiobjAvlItem *>(item); }

/**
 * @brief The biobj logger, which wraps the observed problem.
 *
 * @note Some fields from the observers (coco_observer as well as observer_biobj) are copied here because the
 * observers can be deleted before the logger is finalized and we need these fields for finalization.
 */
class LoggerBiobj final : public Logger {
public:
  /**< @brief The implemented indicators. */

  LoggerBiobj(ObserverBiobj &observer, ProblemPtr inner);

  ~LoggerBiobj() override;

  void EvaluateFunction(std::span<const double> x, std::span<double> y) override;
  void RecommendSolution(std::span<const double> x) override;
  bool FeedSolution(const size_t evaluation, std::span<const double> y);

  /** @brief Saves the information that the algorithm has restarted. */
  void SignalRestart() override;

  /** @brief Returns the indicator with the given index (nullptr if the indicators are not computed). */
  LoggerBiobjIndicator *indicator(const size_t i) const { return indicators_[i].get(); }

private:
  friend class LoggerBiobjIndicator; /* Reads the settings and counters of the logger for its output */

  /** @brief Adds the solution to the archive if it is not dominated, returns whether the archive changed. */
  bool UpdateArchive(LoggerBiobjAvlItem *node_item);

  /** @brief Outputs the solutions in the tree to the given file. Returns the number of solutions. */
  size_t OutputTree(FILE *file, const AvlTree *tree) const;

  /** @brief Outputs the indicator information according to the observer options. */
  void Output(const int update_performed, const LoggerBiobjAvlItem *node_item);

  /** @brief Outputs the final nondominated solutions to the archive file. */
  void WriteFinalArchive();

  /** @brief The decision values are not transformed. */
  bool Inverse(std::span<const double> inner_x, std::span<double> x) const override {
    return IdentityInverse(inner_x, x);
  }

  ObserverBiobjLogNondom log_nondom_mode_{};
  /**< @brief Mode for archiving nondominated solutions. */
  File adat_file_;                /**< @brief File for archiving nondominated solutions (all or final). */
  File mdat_file_;                /**< @brief File for logging recommended solutions */
  bool log_vars_{};               /**< @brief Whether to log the decision values. */
  int precision_x_{};             /**< @brief Precision for outputting decision values. */
  int precision_f_{};             /**< @brief Precision for outputting objective values. */
  bool log_discrete_as_int_{};    /**< @brief Whether to output discrete variables in int or double format. */
  bool algorithm_restarted_{};    /**< @brief Whether the algorithm has restarted (output information to .rdat file). */
  size_t num_func_evaluations_{}; /**< @brief The number of function evaluations performed so far. */
  size_t num_cons_evaluations_{}; /**< @brief The number of evaluations of constraints performed so far. */
  size_t previous_evaluations_{}; /**< @brief The number of evaluations from the previous call to the logger. */
  AvlTreePtr archive_tree_;       /**< @brief The tree keeping currently non-dominated solutions. */
  AvlTreePtr buffer_tree_;        /**< @brief The tree with pointers to nondominated solutions that haven't
                                          been logged yet. */
  /* TODO: Implement other indicators */
  bool compute_indicators_{}; /**< @brief Whether to compute the indicators. */
  std::unique_ptr<LoggerBiobjIndicator> indicators_[LOGGER_BIOBJ_NUMBER_OF_INDICATORS];
};

/**
 * @brief Defines the ordering of AVL tree nodes based on the value of the last objective.
 *
 * @note This ordering is used by the archive_tree.
 */
static int AvlTreeCompareByLastObjective(const void *a, const void *b, void *) {
  const auto *item1 = static_cast<const LoggerBiobjAvlItem *>(a);
  const auto *item2 = static_cast<const LoggerBiobjAvlItem *>(b);
  if (DoubleAlmostEqual(item1->normalized_y[1], item2->normalized_y[1], MO_PRECISION))
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
static int AvlTreeCompareByEvalNumber(const void *a, const void *b, void *) {
  const auto *item1 = static_cast<const LoggerBiobjAvlItem *>(a);
  const auto *item2 = static_cast<const LoggerBiobjAvlItem *>(b);
  if (item1->evaluation_number < item2->evaluation_number)
    return -1;
  else if (item1->evaluation_number > item2->evaluation_number)
    return 1;
  else
    return 0;
}

size_t LoggerBiobj::OutputTree(FILE *file, const AvlTree *tree) const {

  AvlNode *solution;
  size_t i;
  size_t number_of_nodes = 0;

  if (tree->tail) {
    /* There is at least a solution in the tree to output */
    solution = tree->head;
    while (solution != nullptr) {
      fprintf(file, "%zu\t", (static_cast<LoggerBiobjAvlItem *>(solution->item))->evaluation_number);
      for (size_t j = 0; j < number_of_objectives_; j++)
        fprintf(file, "%.*e\t", precision_f_, (static_cast<LoggerBiobjAvlItem *>(solution->item))->y[j]);
      if (log_vars_) {
        for (i = 0; i < number_of_variables_; i++)
          if ((i < number_of_integer_variables_) && (log_discrete_as_int_))
            fprintf(file, "%d\t", DoubleToInt((static_cast<LoggerBiobjAvlItem *>(solution->item))->x[i]));
          else
            fprintf(file, "%.*e\t", precision_x_, (static_cast<LoggerBiobjAvlItem *>(solution->item))->x[i]);
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
 * @return Whether the update was performed.
 */
bool LoggerBiobj::UpdateArchive(LoggerBiobjAvlItem *node_item) {

  AvlNode *node, *next_node, *new_node;
  bool trigger_update = false;
  int dominance;
  bool previous_unavailable = false;

  /* If the node contains an infeasible solution, exit immediately (do not update the tree) */
  if (node_item->is_feasible == 0)
    return false;

  /* Find the first point that is not worse than the new point (nullptr if such point does not exist) */
  node = AvlItemSearchRight(archive_tree_.get(), node_item, nullptr);

  if (node == nullptr) {
    /* The new point is an extreme point */
    trigger_update = true;
    next_node = archive_tree_->head;
  } else {
    dominance = MoGetDominance(node_item->normalized_y, (static_cast<LoggerBiobjAvlItem *>(node->item))->normalized_y);
    if (dominance > -1) {
      trigger_update = true;
      next_node = node->next;
      if (dominance == 1) {
        /* The new point dominates the next point, remove the next point */
        if (compute_indicators_) {
          for (size_t i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
            indicators_[i]->RemoveContribution(
                (static_cast<LoggerBiobjAvlItem *>(node->item))->indicator_contribution[i]);
          }
        }
        AvlItemDelete(buffer_tree_.get(), node->item);
        AvlNodeDelete(archive_tree_.get(), node);
      }
    } else {
      /* The new point is dominated or equal to an existing one, nothing more to do */
      trigger_update = false;
    }
  }

  if (!trigger_update) {
    delete node_item;
  } else {
    /* Perform tree update */
    while (next_node != nullptr) {
      /* Check the dominance relation between the new node and the next node. There are only two possibilities:
       * dominance = 0: the new node and the next node are nondominated
       * dominance = 1: the new node dominates the next node */
      node = next_node;
      dominance =
          MoGetDominance(node_item->normalized_y, (static_cast<LoggerBiobjAvlItem *>(node->item))->normalized_y);
      if (dominance == 1) {
        /* The new point dominates the next point, remove the next point */
        if (compute_indicators_) {
          for (size_t i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
            indicators_[i]->RemoveContribution(
                (static_cast<LoggerBiobjAvlItem *>(node->item))->indicator_contribution[i]);
          }
        }
        next_node = node->next;
        AvlItemDelete(buffer_tree_.get(), node->item);
        AvlNodeDelete(archive_tree_.get(), node);
      } else {
        break;
      }
    }

    new_node = AvlItemInsert(archive_tree_.get(), node_item);
    assert(new_node != nullptr);
    AvlItemInsert(buffer_tree_.get(), node_item);

    if (compute_indicators_) {
      if (node_item->within_ROI) {
        /* Compute indicator value for new node and update the indicator value of the affected nodes */
        LoggerBiobjAvlItem *next_item, *previous_item;

        if (new_node->next != nullptr) {
          next_item = static_cast<LoggerBiobjAvlItem *>(new_node->next->item);
          if (next_item->within_ROI) {
            for (size_t i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
              indicators_[i]->RemoveContribution(next_item->indicator_contribution[i]);
              if (indicators_[i]->name() == "hyp") {
                next_item->indicator_contribution[i] =
                    (node_item->normalized_y[0] - next_item->normalized_y[0]) * (1 - next_item->normalized_y[1]);
                assert(next_item->indicator_contribution[i] >= 0);
              } else {
                coco_error(
                    "logger_biobj_t::UpdateArchive(): Indicator computation not implemented yet for indicator %s",
                    indicators_[i]->name().c_str());
              }
              indicators_[i]->AddContribution(next_item->indicator_contribution[i]);
            }
          }
        }

        previous_unavailable = false;
        if (new_node->prev != nullptr) {
          previous_item = static_cast<LoggerBiobjAvlItem *>(new_node->prev->item);
          if (previous_item->within_ROI) {
            for (size_t i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
              if (indicators_[i]->name() == "hyp") {
                node_item->indicator_contribution[i] =
                    (previous_item->normalized_y[0] - node_item->normalized_y[0]) * (1 - node_item->normalized_y[1]);
                assert(node_item->indicator_contribution[i] >= 0);
              } else {
                coco_error(
                    "logger_biobj_t::UpdateArchive(): Indicator computation not implemented yet for indicator %s",
                    indicators_[i]->name().c_str());
              }
            }
          } else {
            previous_unavailable = true;
          }
        } else {
          previous_unavailable = true;
        }

        if (previous_unavailable) {
          /* Previous item does not exist or is out of ROI, use reference point instead */
          for (size_t i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
            if (indicators_[i]->name() == "hyp") {
              node_item->indicator_contribution[i] =
                  (1 - node_item->normalized_y[0]) * (1 - node_item->normalized_y[1]);
              assert(node_item->indicator_contribution[i] >= 0);
            } else {
              coco_error("logger_biobj_t::UpdateArchive(): Indicator computation not implemented yet for indicator %s",
                         indicators_[i]->name().c_str());
            }
          }
        }

        for (size_t i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
          if (indicators_[i]->name() == "hyp") {
            assert(node_item->indicator_contribution[i] >= 0);
            indicators_[i]->AddContribution(node_item->indicator_contribution[i]);
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
void LoggerBiobjIndicator::InitializeFile(const LoggerBiobj &logger, const ObserverBiobj &observer,
                                          const Problem &problem, File &f, const char *file_ending,
                                          const int output_targets) {
  static const char *header_w_targets = "%%\n"
                                        "%% index = %zu, name = %s\n"
                                        "%% instance = %zu, reference value = %.*e\n"
                                        "%% function evaluation | indicator value | target hit\n";
  static const char *header_wo_targets = "%%\n"
                                         "%% index = %zu, name = %s\n"
                                         "%% instance = %zu, reference value = %.*e\n"
                                         "%% function evaluation | indicator value\n";

  /* Create and open the file */
  const std::string folder_name = JoinPath(observer.result_folder, problem.type());
  CreateDirectory(folder_name);
  const std::string prefix = RemoveFromString(problem.id(), "_i", "_d");
  const std::string path_name =
      JoinPath(folder_name, StringFormat("%s_%s.%s", prefix.c_str(), name_.c_str(), file_ending));
  if (!f.Open(path_name, "a")) {
    coco_error("logger_biobj_indicator_t::InitializeFile() failed to open file '%s'.", path_name.c_str());
  }

  /* Output header */
  if (output_targets)
    fprintf(f, header_w_targets, problem.suite_dep_index(), problem.name().c_str(), problem.suite_dep_instance(),
            logger.precision_f_, best_value_);
  else
    fprintf(f, header_wo_targets, problem.suite_dep_index(), problem.name().c_str(), problem.suite_dep_instance(),
            logger.precision_f_, best_value_);
}

/**
 * @brief Initializes the indicator with name indicator_name.
 *
 * Opens files for writing and resets counters.
 */
LoggerBiobjIndicator::LoggerBiobjIndicator(const LoggerBiobj &logger, const ObserverBiobj &observer,
                                           const Problem &problem, std::string_view name)
    : name_(name) {

  const Suite *suite;
  int info_file_exists = 0;

  coco_debug("Started LoggerBiobjIndicator()");

  assert(problem.suite() != nullptr);
  suite = problem.suite();
  target_hit_ = 0;
  evaluation_logged_ = 0;
  current_value_ = 0;
  additional_penalty_ = DBL_MAX;
  overall_value_ = 0;
  previous_value_ = 0;

  /* The best known value of the indicator */
  if (this->name_ == "hyp")
    SuiteBiobjGetBestHypValue(suite->has_known_optima(), problem.id(), &best_value_);
  else
    coco_error("LoggerBiobjIndicator(): indicator %s not supported", this->name_.c_str());

  targets_ = ObserverTargets(suite->has_known_optima(), observer.lin_target_precision, observer.number_target_triggers,
                             observer.log_target_precision);
  evaluations_ = ObserverEvaluations(observer.number_evaluation_triggers, observer.base_evaluation_triggers,
                                     problem.number_of_variables());

  /* Prepare the info file */
  CreateDirectory(observer.result_folder);
  const std::string path_name =
      JoinPath(observer.result_folder, StringFormat("%s_%s.info", problem.type().c_str(), this->name_.c_str()));
  info_file_exists = FileExists(path_name.c_str());
  info_file_.Open(path_name, "a");
  if (info_file_ == nullptr) {
    coco_error("LoggerBiobjIndicator() failed to open file '%s'.", path_name.c_str());
    return; /* Never reached */
  }

  /* Output header information to the info file */
  const std::string prefix = RemoveFromString(problem.id(), "_i", "_d");
  if (!info_file_exists) {
    /* Output algorithm name */
    assert(problem.suite() != nullptr);
    fprintf(info_file_,
            "suite = '%s', algorithm = '%s', indicator = '%s', folder = '%s', coco_version = '%s', "
            "settings = '%s' \n%% %s",
            problem.suite()->name().c_str(), observer.algorithm_name.c_str(), this->name_.c_str(),
            problem.type().c_str(), coco_version, observer.settings.c_str(), observer.algorithm_info.c_str());
    if (logger.log_nondom_mode_ == LOG_NONDOM_READ)
      fprintf(info_file_, " (reconstructed)");
  }
  if ((observer.previous_function_ != static_cast<long>(problem.suite_dep_function())) ||
      (observer.previous_dimension_ != static_cast<long>(problem.number_of_variables()))) {
    fprintf(info_file_, "\nfunction = %2zu, ", problem.suite_dep_function());
    fprintf(info_file_, "dim = %2zu, ", problem.number_of_variables());
    fprintf(info_file_, "%s_%s.dat", prefix.c_str(), this->name_.c_str());
  }

  /* Initialize the .dat, .tdat and .rdat files */
  InitializeFile(logger, observer, problem, dat_file_, "dat", 1);
  InitializeFile(logger, observer, problem, tdat_file_, "tdat", 0);
  InitializeFile(logger, observer, problem, rdat_file_, "rdat", 0);

  coco_debug("Ended   LoggerBiobjIndicator()");
}

/**
 * @brief Outputs the final information about this indicator.
 */
void LoggerBiobjIndicator::Finalize(const LoggerBiobj &logger) {

  /* Log the last eval_number in the dat file if wasn't already logged */
  if (!target_hit_) {
    fprintf(dat_file_, "%zu\t%.*e\t%.*e\n", logger.num_func_evaluations_, logger.precision_f_, overall_value_,
            logger.precision_f_, targets_.LastTarget());
  }

  /* Log the last eval_number in the tdat file if wasn't already logged */
  if (!evaluation_logged_) {
    fprintf(tdat_file_, "%zu\t%.*e\n", logger.num_func_evaluations_, logger.precision_f_, overall_value_);
  }

  /* Log the information in the info file */
  fprintf(info_file_, ", %zu:%zu|%.1e", logger.suite_dep_instance(), logger.num_func_evaluations_, overall_value_);
  fflush(info_file_);
}

/**
 * @brief Updates and outputs the indicator, see LoggerBiobj::output().
 */
void LoggerBiobjIndicator::Output(const LoggerBiobj &logger, const int update_performed,
                                  const LoggerBiobjAvlItem *node_item) {
  target_hit_ = 0;
  previous_value_ = overall_value_;

  /* If the update was performed, update the overall indicator value */
  if (update_performed) {
    /* Compute the overall_value of an indicator */
    if (name_ == "hyp") {
      if (DoubleAlmostEqual(current_value_, 0, MO_PRECISION)) {
        /* Update the additional penalty for hypervolume (the minimal distance from the nondominated set
         * to the ROI) */
        double new_distance = MoGetDistanceToRoi(node_item->normalized_y);
        additional_penalty_ = DoubleMin(additional_penalty_, new_distance);
        assert(additional_penalty_ >= 0);
      } else {
        additional_penalty_ = 0;
      }
      overall_value_ = best_value_ - current_value_ + additional_penalty_;
    } else {
      coco_error("logger_biobj_indicator_t::output(): Indicator computation not implemented yet for indicator %s",
                 name_.c_str());
    }

    /* Check whether a target was hit */
    target_hit_ = targets_.Trigger(overall_value_);
  } else if ((logger.num_func_evaluations_ == 1) && (node_item->is_feasible == 0)) {
    /* Special case if the solution is infeasible and this is the first evaluation */
    overall_value_ = INFINITY_FOR_LOGGING;
    target_hit_ = targets_.Trigger(overall_value_);
  }

  /* Log to the dat file if a performance target was hit */
  if (target_hit_) {
    fprintf(dat_file_, "%zu\t%.*e\t%.*e\n", logger.num_func_evaluations_, logger.precision_f_, overall_value_,
            logger.precision_f_, targets_.LastTarget());
  }

  if (logger.log_nondom_mode_ == LOG_NONDOM_READ) {
    /* Log to the tdat file the previous indicator value if any evaluation number between the previous and
     * this one matches one of the predefined evaluation numbers. */
    for (size_t j = logger.previous_evaluations_ + 1; j < logger.num_func_evaluations_; j++) {
      evaluation_logged_ = evaluations_.Trigger(j);
      if (evaluation_logged_) {
        fprintf(tdat_file_, "%zu\t%.*e\n", j, logger.precision_f_, previous_value_);
      }
    }
  }

  /* Log to the tdat file if the number of evaluations matches one of the predefined numbers */
  evaluation_logged_ = evaluations_.Trigger(logger.num_func_evaluations_);
  if (evaluation_logged_) {
    fprintf(tdat_file_, "%zu\t%.*e\n", logger.num_func_evaluations_, logger.precision_f_, overall_value_);
  }

  /* Log to the rdat file if the algorithm was restarted */
  if (logger.algorithm_restarted_) {
    fprintf(rdat_file_, "%zu\t%.*e\n", logger.num_func_evaluations_, logger.precision_f_, overall_value_);
  }
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
void LoggerBiobj::Output(const int update_performed, const LoggerBiobjAvlItem *node_item) {

  coco_debug("Started logger_biobj_t::output()");

  if (compute_indicators_) {
    for (size_t i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++)
      indicators_[i]->Output(*this, update_performed, node_item);
  }

  if (algorithm_restarted_)
    algorithm_restarted_ = 0;

  coco_debug("Ended   logger_biobj_t::output()");
}

/**
 * @brief Evaluates the function, increases the number of evaluations and outputs information according to
 * observer options.
 */
void LoggerBiobj::EvaluateFunction(std::span<const double> x, std::span<double> y) {
  LoggerBiobjAvlItem *node_item;

  /* Evaluate the objectives */
  inner_problem_->EvaluateFunction(x, y);
  num_func_evaluations_++;

  /* Evaluate the constraints */
  std::vector<double> constraints(number_of_constraints_);
  if (number_of_constraints_ > 0)
    inner_problem_->EvaluateConstraint(x, constraints, 0);
  num_cons_evaluations_ = evaluations_constraints();

  node_item = new LoggerBiobjAvlItem(*inner_problem_, x, y, constraints, num_func_evaluations_);

  /* Update the archive with the new solution, if it is not dominated by or equal to existing solutions in
   * the archive */
  int update_performed = UpdateArchive(node_item);

  /* If the archive was updated and you need to log all nondominated solutions, output the new solution to
   * nondom_file */
  if (update_performed && (log_nondom_mode_ == LOG_NONDOM_ALL)) {
    OutputTree(adat_file_, buffer_tree_.get());
    AvlTreePurge(buffer_tree_.get());

    /* Flush output so that impatient users can see progress. */
    fflush(adat_file_);
  }

  /* Output according to observer options */
  Output(update_performed, node_item);
}

/**
 * @brief Evaluates the function and outputs information according to observer options to the file with
 * recommendations. The evaluation result is not returned and the evaluation counter is not increased.
 */
void LoggerBiobj::RecommendSolution(std::span<const double> x) {
  size_t i;

  /* Evaluate the objectives */
  std::vector<double> y(number_of_objectives_);
  coco_evaluate_function(inner_problem_.get(), x.data(), y.data()); /* checks x, which recommendations do not */

  /* Evaluate the constraints */
  std::vector<double> constraints(number_of_constraints_);
  if (number_of_constraints_ > 0)
    inner_problem_->EvaluateConstraint(x, constraints, 0);
  num_cons_evaluations_ = evaluations_constraints();

  /* If this is the first call to the recommend function (on any problem), warn about the large
  resulting files (the observer remembers the warning, unless it has already been freed) */
  auto *observer = static_cast<ObserverBiobj *>(this->observer());
  if ((observer != nullptr) && (observer->large_mdat_file_warning_ == 0)) {
    coco_warning("Using the logger_biobj_recommend function may result in the creation of large "
                 "mdat files since every recommended solution is logged");
    observer->large_mdat_file_warning_ = 1;
  }

  /* Log to the mdat file */
  fprintf(mdat_file_, "%zu\t", num_func_evaluations_);
  for (size_t j = 0; j < number_of_objectives_; j++)
    fprintf(mdat_file_, "%.*e\t", precision_f_, y[j]);
  for (size_t j = 0; j < number_of_constraints_; j++)
    fprintf(mdat_file_, "%.*e\t", precision_f_, constraints[j]);
  if (log_vars_) {
    for (i = 0; i < number_of_variables_; i++)
      if ((i < number_of_integer_variables_) && (log_discrete_as_int_))
        fprintf(mdat_file_, "%d\t", DoubleToInt(x[i]));
      else
        fprintf(mdat_file_, "%.*e\t", precision_x_, x[i]);
  }
  fprintf(mdat_file_, "\n");
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
 * @return Whether the archive was updated.
 */
bool LoggerBiobj::FeedSolution(const size_t evaluation, std::span<const double> y) {
  LoggerBiobjAvlItem *node_item;
  bool update_performed;

  assert(log_nondom_mode_ == LOG_NONDOM_READ);

  /* Set the number of evaluations */
  previous_evaluations_ = num_func_evaluations_;
  if (previous_evaluations_ >= evaluation)
    coco_error("logger_biobj_t::FeedSolution(): Evaluation %zu came before evaluation %zu. Note that "
               "the evaluations need to be always increasing.",
               previous_evaluations_, evaluation);
  num_func_evaluations_ = evaluation;

  /* Update the archive with the new solution */
  std::vector<double> x(number_of_variables_);
  for (size_t i = 0; i < number_of_variables_; i++)
    x[i] = 0;
  node_item = new LoggerBiobjAvlItem(*inner_problem_, x, y, {}, num_func_evaluations_);

  /* Update the archive */
  update_performed = UpdateArchive(node_item);

  /* Output according to observer options */
  Output(update_performed, node_item);

  return update_performed;
}

/**
 * @brief Outputs the final nondominated solutions to the archive file.
 */
void LoggerBiobj::WriteFinalArchive() {

  AvlTree *resorted_tree;
  AvlNode *solution;

  coco_debug("Started logger_biobj_t::WriteFinalArchive()");

  /* Re-sort archive_tree according to time stamp and then output it */
  resorted_tree = AvlTreeConstruct(AvlTreeCompareByEvalNumber, nullptr);

  if (archive_tree_->tail) {
    /* There is at least a solution in the tree to output */
    solution = archive_tree_->head;
    while (solution != nullptr) {
      AvlItemInsert(resorted_tree, solution->item);
      solution = solution->next;
    }
  }

  OutputTree(adat_file_, resorted_tree);

  AvlTreeDestruct(resorted_tree);

  coco_debug("Ended   logger_biobj_t::WriteFinalArchive()");
}

/**
 * @brief Outputs the final information and closes the files.
 */
LoggerBiobj::~LoggerBiobj() {

  coco_debug("Started logger_biobj_t::~LoggerBiobj()");

  if (log_nondom_mode_ == LOG_NONDOM_FINAL) {
    WriteFinalArchive();
  }

  if (compute_indicators_) {
    for (size_t i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++) {
      indicators_[i]->Finalize(*this);
      indicators_[i].reset();
    }
  }

  if (((log_nondom_mode_ == LOG_NONDOM_ALL) || (log_nondom_mode_ == LOG_NONDOM_FINAL)) && (adat_file_ != nullptr)) {
    fprintf(adat_file_, "%% evaluations = %zu\n", num_func_evaluations_);
    adat_file_.Close();
  }

  coco_debug("Ended   logger_biobj_t::~LoggerBiobj()");
}

/**
 * @brief Saves the information that the algorithm has restarted
 */
void LoggerBiobj::SignalRestart() {
  if (num_func_evaluations_ > 0)
    algorithm_restarted_ = 1;
}

/**
 * @brief Initializes the biobjective logger.
 *
 * Copies all observer field values that are needed after initialization into logger field values
 * because if the observer is deleted before the suite, the observer is not available anymore when
 * the logger is finalized.
 */
LoggerBiobj::LoggerBiobj(ObserverBiobj &observer, ProblemPtr inner) : Logger(observer, std::move(inner)) {
  const char nondom_folder_name[] = "archive";

  coco_debug("Started LoggerBiobj()");

  if (inner_problem_->number_of_objectives() != 2) {
    coco_error("LoggerBiobj(): The bi-objective logger cannot log a problem with %zu objective(s)",
               inner_problem_->number_of_objectives());
    return; /* Never reached */
  }


  num_func_evaluations_ = 0;
  num_cons_evaluations_ = 0;
  previous_evaluations_ = 0;
  algorithm_restarted_ = 0;

  /* Copy values from the observer that you might need even if it does not exist any more */
  log_nondom_mode_ = observer.log_nondom_mode_;
  compute_indicators_ = observer.compute_indicators_;
  precision_x_ = observer.precision_x;
  precision_f_ = observer.precision_f;
  log_discrete_as_int_ = observer.log_discrete_as_int;

  if (((observer.log_vars_mode_ == LOG_VARS_LOW_DIM) && (inner_problem_->number_of_variables() > 5)) ||
      (observer.log_vars_mode_ == LOG_VARS_NEVER))
    log_vars_ = 0;
  else
    log_vars_ = 1;

  /* Initialize logging of nondominated solutions into the archive file */
  if ((log_nondom_mode_ == LOG_NONDOM_ALL) || (log_nondom_mode_ == LOG_NONDOM_FINAL)) {

    /* Create the path to the file */
    const std::string folder_name = JoinPath(observer.result_folder, nondom_folder_name);
    CreateDirectory(folder_name);

    /* Construct file name */
    const std::string path_name =
        JoinPath(folder_name, inner_problem_->id() +
                                  ((log_nondom_mode_ == LOG_NONDOM_ALL) ? "_nondom_all.adat" : "_nondom_final.adat"));

    /* Open and initialize the archive file */
    adat_file_.Open(path_name, "a");
    if (adat_file_ == nullptr) {
      coco_error("LoggerBiobj() failed to open file '%s'.", path_name.c_str());
      return; /* Never reached */
    }

    /* Output header information */
    fprintf(adat_file_, "%% instance = %zu, name = %s\n", inner_problem_->suite_dep_instance(),
            inner_problem_->name().c_str());
    if (log_vars_) {
      fprintf(adat_file_, "%% function evaluation | %zu objectives | %zu variables\n",
              inner_problem_->number_of_objectives(), inner_problem_->number_of_variables());
    } else {
      fprintf(adat_file_, "%% function evaluation | %zu objectives \n", inner_problem_->number_of_objectives());
    }
  }

  /* Initialize the AVL trees */
  archive_tree_.reset(AvlTreeConstruct(AvlTreeCompareByLastObjective, LoggerBiobjNodeFree));
  buffer_tree_.reset(AvlTreeConstruct(AvlTreeCompareByEvalNumber, nullptr));

  /* Initialize the indicators */
  if (compute_indicators_) {
    for (size_t i = 0; i < LOGGER_BIOBJ_NUMBER_OF_INDICATORS; i++)
      indicators_[i] =
          std::make_unique<LoggerBiobjIndicator>(*this, observer, *inner_problem_, LOGGER_BIOBJ_INDICATORS[i]);

    observer.previous_function_ = static_cast<long>(inner_problem_->suite_dep_function());
    observer.previous_dimension_ = static_cast<long>(inner_problem_->number_of_variables());
  }

  /* Initialize the file with recommendations */

  /* Create the path to the file */
  const std::string folder_name = JoinPath(observer.result_folder, type());
  CreateDirectory(folder_name);

  /* Construct file name */
  const std::string prefix = RemoveFromString(id(), "_i", "_d");
  const std::string path_name = JoinPath(folder_name, prefix + ".mdat");

  /* Open and initialize the recommendation file */
  mdat_file_.Open(path_name, "a");
  if (mdat_file_ == nullptr) {
    coco_error("LoggerBiobj() failed to open file '%s'.", path_name.c_str());
    return; /* Never reached */
  }

  /* Output header information */
  fprintf(mdat_file_, "%% instance = %zu, name = %s\n%% function evaluation | %zu objectives",
          inner_problem_->suite_dep_instance(), inner_problem_->name().c_str(), inner_problem_->number_of_objectives());
  if (inner_problem_->number_of_constraints() > 0)
    fprintf(mdat_file_, " | %zu constraints", inner_problem_->number_of_constraints());
  if (log_vars_)
    fprintf(mdat_file_, " | %zu variables", inner_problem_->number_of_variables());
  fprintf(mdat_file_, "\n");

  coco_debug("Ended   LoggerBiobj()");
}

ProblemPtr ObserverBiobj::AddLogger(ProblemPtr problem) {
  if (number_of_loggers() > 0)
    coco_error("LoggerBiobj(): The observed problem must be closed before a new problem can be observed");
  return std::make_unique<LoggerBiobj>(*this, std::move(problem));
}

} // namespace coco

/* C interface (declared in coco.h) */

/**
 * @brief Feeds the solution to the bi-objective logger of the given problem, see
 * LoggerBiobj::FeedSolution().
 */
int coco_logger_biobj_feed_solution(coco_problem_t *problem, const size_t evaluation, const double *y) {
  auto *logger = dynamic_cast<coco::LoggerBiobj *>(FromC(problem));
  if (logger == nullptr)
    coco_error("coco_logger_biobj_feed_solution(): Problem %s is not observed by a bbob-biobj logger",
               problem == nullptr ? "NULL" : FromC(problem)->id().c_str());
  return logger->FeedSolution(evaluation, {y, logger->number_of_objectives()});
}
