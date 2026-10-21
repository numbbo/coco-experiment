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

/**
 * @brief The COCO archive structure.
 *
 * The archive structure is used for pre-processing archives of non-dominated solutions.
 */
struct coco_archive_s {
public:
  /** @brief Creates the archive of the given problem, holding its two extreme points. */
  coco_archive_s(std::string_view suite_name, const size_t function, const size_t dimension, const size_t instance);
  ~coco_archive_s();

  coco_archive_s(const coco_archive_s &) = delete;
  coco_archive_s &operator=(const coco_archive_s &) = delete;

  /** @brief Adds the solution to the archive if it is not dominated, returns whether the archive was updated. */
  bool add_solution(const double y1, const double y2, std::string_view text);

  /** @brief Returns the text of the next solution (the extremes first) or "" after the last one. */
  const char *get_next_solution_text();

  size_t get_number_of_solutions();
  double get_hypervolume();

private:
  /** @brief Updates the fields returned by the getters. */
  void update();

  avl_tree_t *tree;               /**< @brief The AVL tree with non-dominated solutions. */
  std::vector<double> ideal;      /**< @brief The ideal point. */
  std::vector<double> nadir;      /**< @brief The nadir point. */
  bool is_up_to_date;             /**< @brief Whether archive fields have been updated since last addition. */
  size_t number_of_solutions;     /**< @brief Number of solutions in the archive. */
  double hypervolume;             /**< @brief Hypervolume of the solutions in the archive. */
  avl_node_t *current_solution;   /**< @brief Current solution (to return). */
  avl_node_t *extreme1;           /**< @brief Pointer to the first extreme solution. */
  avl_node_t *extreme2;           /**< @brief Pointer to the second extreme solution. */
  bool extremes_already_returned; /**< @brief Whether the extreme solutions have already been returned. */
};

/**
 * @brief The type for the node's item in the AVL tree used by the archive.
 *
 * Contains information on the rounded normalized objective values (normalized_y), which are used for
 * computing the indicators and the text, which is used for output.
 */
struct coco_archive_avl_item_t {
  coco_archive_avl_item_t(std::span<const double> y, std::span<const double> ideal, std::span<const double> nadir,
                          std::string_view text)
      : normalized_y(mo_normalize(y, ideal, nadir)), text(text) {}

  std::vector<double> normalized_y; /**< @brief The values of normalized objectives of this solution. */
  std::string text;                 /**< @brief The text describing the solution (the whole line of the archive). */
};

/**
 * @brief Creates and returns the information on the solution in the form of a node's item in the AVL tree.
 */
/**
 * @brief Frees the data of the given coco_archive_avl_item_t.
 */
static void coco_archive_node_item_free(void *item, void *) { delete static_cast<coco_archive_avl_item_t *>(item); }

/**
 * @brief Defines the ordering of AVL tree nodes based on the value of the last objective.
 */
static int coco_archive_compare_by_last_objective(const void *a, const void *b, void *) {
  const auto *item1 = static_cast<const coco_archive_avl_item_t *>(a);
  const auto *item2 = static_cast<const coco_archive_avl_item_t *>(b);
  if (coco_double_almost_equal(item1->normalized_y[1], item2->normalized_y[1], mo_precision))
    return 0;
  else if (item1->normalized_y[1] < item2->normalized_y[1])
    return -1;
  else
    return 1;
}

coco_archive_s::coco_archive_s(std::string_view suite_name, const size_t function, const size_t dimension,
                               const size_t instance)
    : is_up_to_date(0), number_of_solutions(0), hypervolume(0.0), current_solution(nullptr), extreme1(nullptr),
      extreme2(nullptr), extremes_already_returned(0) {

  int output_precision = 15;
  coco_suite_t *suite;
  const std::string suite_instance = coco_string_format("instances: %zu", instance);
  const std::string suite_options = coco_string_format("dimensions: %zu function_indices: %zu", dimension, function);
  coco_problem_t *problem;
  std::string text;

  /* Initialize the AVL tree */
  tree = avl_tree_construct(coco_archive_compare_by_last_objective, coco_archive_node_item_free);

  suite = coco_suite_create(suite_name, suite_instance, suite_options);
  if (suite == nullptr) {
    coco_error("coco_archive(): cannot create suite '%.*s'", static_cast<int>(suite_name.size()), suite_name.data());
    return; /* Never reached */
  }
  problem = coco_suite_get_next_problem(suite, nullptr);
  if (problem == nullptr) {
    coco_error("coco_archive(): cannot create problem f%02zu_i%02zu_d%02zu in suite '%.*s'", function, instance,
               dimension, static_cast<int>(suite_name.size()), suite_name.data());
    return; /* Never reached */
  }

  /* Store the ideal and nadir points */
  ideal.assign(problem->get_best_value().begin(), problem->get_best_value().begin() + 2);
  nadir.assign(problem->get_nadir_value().begin(), problem->get_nadir_value().begin() + 2);

  /* Add the extreme points to the archive */
  text = coco_string_format("0\t%.*e\t%.*e\n", output_precision, nadir[0], output_precision, ideal[1]);
  [[maybe_unused]] int updated = add_solution(nadir[0], ideal[1], text);
  assert(updated == 1);

  text = coco_string_format("0\t%.*e\t%.*e\n", output_precision, ideal[0], output_precision, nadir[1]);
  updated = add_solution(ideal[0], nadir[1], text);
  assert(updated == 1);

  extreme1 = tree->head;
  extreme2 = tree->tail;
  assert(extreme1 != extreme2);

  coco_suite_free(suite);
}

coco_archive_s::~coco_archive_s() { avl_tree_destruct(tree); }

/**
 * The archive always contains the two extreme solutions
 */
coco_archive_t *coco_archive(const char *suite_name, const size_t function, const size_t dimension,
                             const size_t instance) {
  return new coco_archive_t(coco_c_string_view(suite_name), function, dimension, instance);
}

bool coco_archive_s::add_solution(const double y1, const double y2, std::string_view text) {

  coco_archive_avl_item_t *insert_item;
  avl_node_t *node, *next_node;
  bool updated = false;
  int dominance;

  const double y[2] = {y1, y2};
  insert_item = new coco_archive_avl_item_t(y, ideal, nadir, text);

  /* Find the first point that is not worse than the new point (nullptr if such point does not exist) */
  node = avl_item_search_right(tree, insert_item, nullptr);

  if (node == nullptr) {
    /* The new point is an extreme point */
    updated = true;
    next_node = tree->head;
  } else {
    dominance =
        mo_get_dominance(insert_item->normalized_y, (static_cast<coco_archive_avl_item_t *>(node->item))->normalized_y);
    if (dominance > -1) {
      updated = true;
      next_node = node->next;
      if (dominance == 1) {
        /* The new point dominates the next point, remove the next point */
        assert((node != extreme1) && (node != extreme2));
        avl_node_delete(tree, node);
      }
    } else {
      /* The new point is dominated or equal to an existing one, ignore */
      updated = false;
    }
  }

  if (!updated) {
    coco_archive_node_item_free(insert_item, nullptr);
  } else {
    /* Perform tree update */
    while (next_node != nullptr) {
      /* Check the dominance relation between the new node and the next node. There are only two possibilities:
       * dominance = 0: the new node and the next node are nondominated
       * dominance = 1: the new node dominates the next node */
      node = next_node;
      dominance = mo_get_dominance(insert_item->normalized_y,
                                   (static_cast<coco_archive_avl_item_t *>(node->item))->normalized_y);
      if (dominance == 1) {
        next_node = node->next;
        /* The new point dominates the next point, remove the next point */
        assert((node != extreme1) && (node != extreme2));
        avl_node_delete(tree, node);
      } else {
        break;
      }
    }

    if (avl_item_insert(tree, insert_item) == nullptr) {
      coco_archive_node_item_free(insert_item, nullptr);
      coco_warning("Solution %.*s did not update the archive", static_cast<int>(text.size()), text.data());
      updated = false;
    }

    is_up_to_date = 0;
  }

  return updated;
}

int coco_archive_add_solution(coco_archive_t *archive, const double y1, const double y2, const char *text) {
  return archive->add_solution(y1, y2, coco_c_string_view(text));
}

/**
 * @brief Updates the archive fields returned by the getters.
 */
void coco_archive_s::update() {

  double hyp;

  if (!is_up_to_date) {

    avl_node_t *node, *left_node;
    coco_archive_avl_item_t *node_item, *left_node_item;

    /* Updates number_of_solutions */

    number_of_solutions = avl_count(tree);

    /* Updates hypervolume */

    node = tree->head;
    hypervolume = 0; /* Hypervolume of the extreme point equals 0 */
    while (node->next) {
      /* Add hypervolume contributions of the other points that are within ROI */
      left_node = node->next;
      node_item = static_cast<coco_archive_avl_item_t *>(node->item);
      left_node_item = static_cast<coco_archive_avl_item_t *>(left_node->item);
      if (mo_is_within_ROI(left_node_item->normalized_y)) {
        if (mo_is_within_ROI(node_item->normalized_y))
          hyp = (node_item->normalized_y[0] - left_node_item->normalized_y[0]) * (1 - left_node_item->normalized_y[1]);
        else
          hyp = (1 - left_node_item->normalized_y[0]) * (1 - left_node_item->normalized_y[1]);
        assert(hyp >= 0);
        hypervolume += hyp;
      }
      node = left_node;
    }

    is_up_to_date = 1;
    current_solution = nullptr;
    extremes_already_returned = 0;
  }
}

const char *coco_archive_s::get_next_solution_text() {

  const char *text;

  update();

  if (!extremes_already_returned) {

    if (current_solution == nullptr) {
      /* Return the first extreme */
      text = (static_cast<coco_archive_avl_item_t *>(extreme1->item))->text.c_str();
      current_solution = extreme2;
      return text;
    }

    if (current_solution == extreme2) {
      /* Return the second extreme */
      text = (static_cast<coco_archive_avl_item_t *>(extreme2->item))->text.c_str();
      extremes_already_returned = 1;
      current_solution = tree->head;
      return text;
    }

  } else {

    if (current_solution == nullptr)
      return "";

    if ((current_solution == extreme1) || (current_solution == extreme2)) {
      /* Skip this one */
      current_solution = current_solution->next;
      return get_next_solution_text();
    }

    /* Return the current solution and move to the next */
    text = (static_cast<coco_archive_avl_item_t *>(current_solution->item))->text.c_str();
    current_solution = current_solution->next;
    return text;
  }

  return nullptr; /* This point should never be reached. */
}

const char *coco_archive_get_next_solution_text(coco_archive_t *archive) { return archive->get_next_solution_text(); }

size_t coco_archive_s::get_number_of_solutions() {
  update();
  return number_of_solutions;
}

double coco_archive_s::get_hypervolume() {
  update();
  return hypervolume;
}

size_t coco_archive_get_number_of_solutions(coco_archive_t *archive) { return archive->get_number_of_solutions(); }

double coco_archive_get_hypervolume(coco_archive_t *archive) { return archive->get_hypervolume(); }

void coco_archive_free(coco_archive_t *archive) {

  assert(archive != nullptr);
  delete archive;
}
