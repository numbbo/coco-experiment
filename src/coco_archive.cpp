/**
 * @file coco_archive.cpp
 * @brief Definitions of functions regarding COCO archives.
 *
 * COCO archives are used to do some pre-processing on the bi-objective archive files. Namely, through a
 * wrapper written in Python, these functions are used to merge archives and compute their hypervolumes.
 */

#include "coco.h"
#include "coco_utilities.cpp"
#include "mo_utilities.cpp"
#include "mo_avl_tree.cpp"

/** @brief The type declared in coco.h, an empty base class of Archive (see coco_internal.h). */
struct coco_archive_s {};

namespace coco {

/**
 * @brief The COCO archive structure.
 *
 * The archive structure is used for pre-processing archives of non-dominated solutions.
 */
struct Archive : public coco_archive_s {
public:
  /** @brief Creates the archive of the given problem, holding its two extreme points. */
  Archive(std::string_view suite_name, const size_t function, const size_t dimension, const size_t instance);
  ~Archive();

  Archive(const Archive &) = delete;
  Archive &operator=(const Archive &) = delete;

  /** @brief Adds the solution to the archive if it is not dominated, returns whether the archive was updated. */
  bool AddSolution(const double y1, const double y2, std::string_view text);

  /** @brief Returns the text of the next solution (the extremes first) or "" after the last one. */
  const char *GetNextSolutionText();

  size_t GetNumberOfSolutions();
  double GetHypervolume();

private:
  /** @brief Updates the fields returned by the getters. */
  void Update();

  AvlTree *tree_;                  /**< @brief The AVL tree with non-dominated solutions. */
  std::vector<double> ideal_;      /**< @brief The ideal point. */
  std::vector<double> nadir_;      /**< @brief The nadir point. */
  bool is_up_to_date_;             /**< @brief Whether archive fields have been updated since last addition. */
  size_t number_of_solutions_;     /**< @brief Number of solutions in the archive. */
  double hypervolume_;             /**< @brief Hypervolume of the solutions in the archive. */
  AvlNode *current_solution_;      /**< @brief Current solution (to return). */
  AvlNode *extreme1_;              /**< @brief Pointer to the first extreme solution. */
  AvlNode *extreme2_;              /**< @brief Pointer to the second extreme solution. */
  bool extremes_already_returned_; /**< @brief Whether the extreme solutions have already been returned. */
};

/**
 * @brief The type for the node's item in the AVL tree used by the archive.
 *
 * Contains information on the rounded normalized objective values (normalized_y), which are used for
 * computing the indicators and the text, which is used for output.
 */
struct ArchiveAvlItem {
  ArchiveAvlItem(std::span<const double> y, std::span<const double> ideal, std::span<const double> nadir,
                 std::string_view text)
      : normalized_y(MoNormalize(y, ideal, nadir)), text(text) {}

  std::vector<double> normalized_y; /**< @brief The values of normalized objectives of this solution. */
  std::string text;                 /**< @brief The text describing the solution (the whole line of the archive). */
};

/**
 * @brief Creates and returns the information on the solution in the form of a node's item in the AVL tree.
 */
/**
 * @brief Frees the data of the given ArchiveAvlItem.
 */
static void ArchiveNodeItemFree(void *item, void *) { delete static_cast<ArchiveAvlItem *>(item); }

/**
 * @brief Defines the ordering of AVL tree nodes based on the value of the last objective.
 */
static int ArchiveCompareByLastObjective(const void *a, const void *b, void *) {
  const auto *item1 = static_cast<const ArchiveAvlItem *>(a);
  const auto *item2 = static_cast<const ArchiveAvlItem *>(b);
  if (DoubleAlmostEqual(item1->normalized_y[1], item2->normalized_y[1], MO_PRECISION))
    return 0;
  else if (item1->normalized_y[1] < item2->normalized_y[1])
    return -1;
  else
    return 1;
}

Archive::Archive(std::string_view suite_name, const size_t function, const size_t dimension, const size_t instance)
    : is_up_to_date_(0), number_of_solutions_(0), hypervolume_(0.0), current_solution_(nullptr), extreme1_(nullptr),
      extreme2_(nullptr), extremes_already_returned_(0) {

  int output_precision = 15;
  Suite *suite;
  const std::string suite_instance = StringFormat("instances: %zu", instance);
  const std::string suite_options = StringFormat("dimensions: %zu function_indices: %zu", dimension, function);
  Problem *problem;
  std::string text;

  /* Initialize the AVL tree */
  tree_ = AvlTreeConstruct(ArchiveCompareByLastObjective, ArchiveNodeItemFree);

  suite = Suite::Create(suite_name, suite_instance, suite_options);
  if (suite == nullptr) {
    coco_error("coco_archive(): cannot create suite '%.*s'", static_cast<int>(suite_name.size()), suite_name.data());
    return; /* Never reached */
  }
  problem = FromC(coco_suite_get_next_problem(suite, nullptr));
  if (problem == nullptr) {
    coco_error("coco_archive(): cannot create problem f%02zu_i%02zu_d%02zu in suite '%.*s'", function, instance,
               dimension, static_cast<int>(suite_name.size()), suite_name.data());
    return; /* Never reached */
  }

  /* Store the ideal and nadir points */
  ideal_.assign(problem->best_value().begin(), problem->best_value().begin() + 2);
  nadir_.assign(problem->nadir_value().begin(), problem->nadir_value().begin() + 2);

  /* Add the extreme points to the archive */
  text = StringFormat("0\t%.*e\t%.*e\n", output_precision, nadir_[0], output_precision, ideal_[1]);
  [[maybe_unused]] int updated = AddSolution(nadir_[0], ideal_[1], text);
  assert(updated == 1);

  text = StringFormat("0\t%.*e\t%.*e\n", output_precision, ideal_[0], output_precision, nadir_[1]);
  updated = AddSolution(ideal_[0], nadir_[1], text);
  assert(updated == 1);

  extreme1_ = tree_->head;
  extreme2_ = tree_->tail;
  assert(extreme1_ != extreme2_);

  coco_suite_free(suite);
}

Archive::~Archive() { AvlTreeDestruct(tree_); }

bool Archive::AddSolution(const double y1, const double y2, std::string_view text) {

  ArchiveAvlItem *insert_item;
  AvlNode *node, *next_node;
  bool updated = false;
  int dominance;

  const double y[2] = {y1, y2};
  insert_item = new ArchiveAvlItem(y, ideal_, nadir_, text);

  /* Find the first point that is not worse than the new point (nullptr if such point does not exist) */
  node = AvlItemSearchRight(tree_, insert_item, nullptr);

  if (node == nullptr) {
    /* The new point is an extreme point */
    updated = true;
    next_node = tree_->head;
  } else {
    dominance = MoGetDominance(insert_item->normalized_y, (static_cast<ArchiveAvlItem *>(node->item))->normalized_y);
    if (dominance > -1) {
      updated = true;
      next_node = node->next;
      if (dominance == 1) {
        /* The new point dominates the next point, remove the next point */
        assert((node != extreme1_) && (node != extreme2_));
        AvlNodeDelete(tree_, node);
      }
    } else {
      /* The new point is dominated or equal to an existing one, ignore */
      updated = false;
    }
  }

  if (!updated) {
    ArchiveNodeItemFree(insert_item, nullptr);
  } else {
    /* Perform tree update */
    while (next_node != nullptr) {
      /* Check the dominance relation between the new node and the next node. There are only two possibilities:
       * dominance = 0: the new node and the next node are nondominated
       * dominance = 1: the new node dominates the next node */
      node = next_node;
      dominance = MoGetDominance(insert_item->normalized_y, (static_cast<ArchiveAvlItem *>(node->item))->normalized_y);
      if (dominance == 1) {
        next_node = node->next;
        /* The new point dominates the next point, remove the next point */
        assert((node != extreme1_) && (node != extreme2_));
        AvlNodeDelete(tree_, node);
      } else {
        break;
      }
    }

    if (AvlItemInsert(tree_, insert_item) == nullptr) {
      ArchiveNodeItemFree(insert_item, nullptr);
      coco_warning("Solution %.*s did not update the archive", static_cast<int>(text.size()), text.data());
      updated = false;
    }

    is_up_to_date_ = 0;
  }

  return updated;
}

/**
 * @brief Updates the archive fields returned by the getters.
 */
void Archive::Update() {

  double hyp;

  if (!is_up_to_date_) {

    AvlNode *node, *left_node;
    ArchiveAvlItem *node_item, *left_node_item;

    /* Updates number_of_solutions */

    number_of_solutions_ = AvlCount(tree_);

    /* Updates hypervolume */

    node = tree_->head;
    hypervolume_ = 0; /* Hypervolume of the extreme point equals 0 */
    while (node->next) {
      /* Add hypervolume contributions of the other points that are within ROI */
      left_node = node->next;
      node_item = static_cast<ArchiveAvlItem *>(node->item);
      left_node_item = static_cast<ArchiveAvlItem *>(left_node->item);
      if (MoIsWithinRoi(left_node_item->normalized_y)) {
        if (MoIsWithinRoi(node_item->normalized_y))
          hyp = (node_item->normalized_y[0] - left_node_item->normalized_y[0]) * (1 - left_node_item->normalized_y[1]);
        else
          hyp = (1 - left_node_item->normalized_y[0]) * (1 - left_node_item->normalized_y[1]);
        assert(hyp >= 0);
        hypervolume_ += hyp;
      }
      node = left_node;
    }

    is_up_to_date_ = 1;
    current_solution_ = nullptr;
    extremes_already_returned_ = 0;
  }
}

const char *Archive::GetNextSolutionText() {

  const char *text;

  Update();

  if (!extremes_already_returned_) {

    if (current_solution_ == nullptr) {
      /* Return the first extreme */
      text = (static_cast<ArchiveAvlItem *>(extreme1_->item))->text.c_str();
      current_solution_ = extreme2_;
      return text;
    }

    if (current_solution_ == extreme2_) {
      /* Return the second extreme */
      text = (static_cast<ArchiveAvlItem *>(extreme2_->item))->text.c_str();
      extremes_already_returned_ = 1;
      current_solution_ = tree_->head;
      return text;
    }

  } else {

    if (current_solution_ == nullptr)
      return "";

    if ((current_solution_ == extreme1_) || (current_solution_ == extreme2_)) {
      /* Skip this one */
      current_solution_ = current_solution_->next;
      return GetNextSolutionText();
    }

    /* Return the current solution and move to the next */
    text = (static_cast<ArchiveAvlItem *>(current_solution_->item))->text.c_str();
    current_solution_ = current_solution_->next;
    return text;
  }

  return nullptr; /* This point should never be reached. */
}

size_t Archive::GetNumberOfSolutions() {
  Update();
  return number_of_solutions_;
}

double Archive::GetHypervolume() {
  Update();
  return hypervolume_;
}

} // namespace coco

/* C interface (declared in coco.h) */

/** @brief Converts a pointer of the C interface to the internal class. */
inline coco::Archive *FromC(coco_archive_t *archive) { return static_cast<coco::Archive *>(archive); }

/**
 * The archive always contains the two extreme solutions
 */
coco_archive_t *coco_archive(const char *suite_name, const size_t function, const size_t dimension,
                             const size_t instance) {
  return new coco::Archive(coco::CStringView(suite_name), function, dimension, instance);
}

int coco_archive_add_solution(coco_archive_t *archive, const double y1, const double y2, const char *text) {
  return FromC(archive)->AddSolution(y1, y2, coco::CStringView(text));
}

const char *coco_archive_get_next_solution_text(coco_archive_t *archive) {
  return FromC(archive)->GetNextSolutionText();
}

size_t coco_archive_get_number_of_solutions(coco_archive_t *archive) { return FromC(archive)->GetNumberOfSolutions(); }

double coco_archive_get_hypervolume(coco_archive_t *archive) { return FromC(archive)->GetHypervolume(); }

void coco_archive_free(coco_archive_t *archive) {

  assert(archive != nullptr);
  delete FromC(archive);
}
