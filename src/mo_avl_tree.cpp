/*****************************************************************************

 avl.c - Source code for libavl

 Copyright (c) 1998  Michael H. Buselli <cosine@cosine.org>
 Copyright (c) 2000-2009  Wessel Dankers <wsl@fruit.je>

 This file is part of libavl.

 libavl is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation, either version 3 of
 the License, or (at your option) any later version.

 libavl is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU General Public License
 and a copy of the GNU Lesser General Public License along with
 libavl.  If not, see <http://www.gnu.org/licenses/>.

 Augmented AVL-tree. Original by Michael H. Buselli <cosine@cosine.org>.

 Modified by Wessel Dankers <wsl@fruit.je> to add a bunch of bloat
 to the source code, change the interface and replace a few bugs.
 Mail him if you find any new bugs.

 Renamed and additionally modified by BOBBies to fit the COCO platform.

 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

namespace coco {

/* In order to easily (un)comment unused functions */
#define AVL_TREE_COMMENT_UNUSED 1

/* We need either depths, counts or both (the latter being the default) */
#if !defined(AVL_DEPTH) && !defined(AVL_COUNT)
#define AVL_DEPTH
#define AVL_COUNT
#endif

/* User supplied function to compare two items like strcmp() does.
 * For example: compare(a,b) will return:
 *   -1  if a < b
 *    0  if a = b
 *    1  if a > b
 */
using AvlCompare = int (*)(const void *a, const void *b, void *userdata);

/* User supplied function to delete an item when a node is free()d.
 * If nullptr, the item is not free()d.
 */
using AvlFree = void (*)(void *item, void *userdata);

#define AVL_CMP(a, b) ((a) < (b) ? -1 : (a) != (b))

#if defined(AVL_COUNT) && defined(AVL_DEPTH)
#define AVL_NODE_INITIALIZER(item) {0, 0, 0, 0, 0, (item), 0, 0}
#else
#define AVL_NODE_INITIALIZER(item) {0, 0, 0, 0, 0, (item), 0}
#endif

struct AvlNode {
  AvlNode *prev;
  AvlNode *next;
  AvlNode *parent;
  AvlNode *left;
  AvlNode *right;
  void *item;
#ifdef AVL_COUNT
  unsigned long count;
#endif
#ifdef AVL_DEPTH
  unsigned char depth;
#endif
};

struct AvlTree {
  AvlNode *top;
  AvlNode *head;
  AvlNode *tail;
  AvlCompare cmpitem;
  AvlFree freeitem;
  void *userdata;
};

static void AvlRebalance(AvlTree *, AvlNode *);
static AvlNode *AvlNodeInsertAfter(AvlTree *avltree, AvlNode *node, AvlNode *newnode);

#ifdef AVL_COUNT
#define AVL_NODE_COUNT(n) ((n) ? (n)->count : 0)
#define AVL_L_COUNT(n) (AVL_NODE_COUNT((n)->left))
#define AVL_R_COUNT(n) (AVL_NODE_COUNT((n)->right))
#define AVL_CALC_COUNT(n) (AVL_L_COUNT(n) + AVL_R_COUNT(n) + 1)
#endif

#ifdef AVL_DEPTH
#define AVL_NODE_DEPTH(n) ((n) ? (n)->depth : 0)
#define AVL_L_DEPTH(n) (AVL_NODE_DEPTH((n)->left))
#define AVL_R_DEPTH(n) (AVL_NODE_DEPTH((n)->right))
#define AVL_CALC_DEPTH(n)                                                                                              \
  (static_cast<unsigned char>((AVL_L_DEPTH(n) > AVL_R_DEPTH(n) ? AVL_L_DEPTH(n) : AVL_R_DEPTH(n)) + 1))
#endif

static int AvlCheckBalance(AvlNode *avlnode) {
#ifdef AVL_DEPTH
  int d = AVL_R_DEPTH(avlnode) - AVL_L_DEPTH(avlnode);
  return d < -1 ? -1 : d > 1;
#else
  /*  int d;
   *  d = ffs(AVL_R_COUNT(avl_node)) - ffs(AVL_L_COUNT(avl_node));
   *  d = d < -1 ? -1 : d > 1;
   */
#ifdef AVL_COUNT
  int pl, r;

  pl = ffs(AVL_L_COUNT(avlnode));
  r = AVL_R_COUNT(avlnode);

  if (r >> pl + 1)
    return 1;
  if (pl < 2 || r >> pl - 2)
    return 0;
  return -1;
#else
#error No balancing possible.
#endif
#endif
}

#ifdef AVL_COUNT
static unsigned long AvlCount(const AvlTree *avltree) {
  if (!avltree)
    return 0;
  return AVL_NODE_COUNT(avltree->top);
}

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_at(const avl_tree_t *avltree, unsigned long index) {
  avl_node_t *avlnode;
  unsigned long c;

  if (!avltree)
    return nullptr;

  avlnode = avltree->top;

  while (avlnode) {
    c = AVL_L_COUNT(avlnode);

    if (index < c) {
      avlnode = avlnode->left;
    } else if (index > c) {
      avlnode = avlnode->right;
      index -= c + 1;
    } else {
      return avlnode;
    }
  }
  return nullptr;
}
#endif

#if (!AVL_TREE_COMMENT_UNUSED)
static unsigned long avl_index(const avl_node_t *avlnode) {
  avl_node_t *next;

  if (!avlnode)
    return 0;

  unsigned long c = AVL_L_COUNT(avlnode);

  while ((next = avlnode->parent)) {
    if (avlnode == next->right)
      c += AVL_L_COUNT(next) + 1;
    avlnode = next;
  }

  return c;
}
#endif
#endif

#if (!AVL_TREE_COMMENT_UNUSED)
static const avl_node_t *avl_search_leftmost_equal(const avl_tree_t *tree, const avl_node_t *node, const void *item) {
  avl_compare_t cmp = tree->cmpitem;
  void *userdata = tree->userdata;
  const avl_node_t *r = node;

  for (;;) {
    for (;;) {
      node = node->left;
      if (!node)
        return r;
      if (cmp(item, node->item, userdata))
        break;
      r = node;
    }
    for (;;) {
      node = node->right;
      if (!node)
        return r;
      if (!cmp(item, node->item, userdata))
        break;
    }
    r = node;
  }

  return nullptr; /* To silence the compiler */
}
#endif

static const AvlNode *AvlSearchRightmostEqual(const AvlTree *tree, const AvlNode *node, const void *item) {
  AvlCompare cmp = tree->cmpitem;
  void *userdata = tree->userdata;
  const AvlNode *r = node;

  for (;;) {
    for (;;) {
      node = node->right;
      if (!node)
        return r;
      if (cmp(item, node->item, userdata))
        break;
      r = node;
    }
    for (;;) {
      node = node->left;
      if (!node)
        return r;
      if (!cmp(item, node->item, userdata))
        break;
    }
    r = node;
  }

  return nullptr; /* To silence the compiler */
}

/* Searches for an item, returning either some exact
 * match, or (if no exact match could be found) the first (leftmost)
 * of the nodes that have an item larger than the search item.
 * If exact is not nullptr, *exact will be set to:
 *    0  if the returned node is unequal or nullptr
 *    1  if the returned node is equal
 * Returns nullptr if no equal or larger element could be found.
 * O(lg n) */
#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_search_leftish(const avl_tree_t *tree, const void *item, int *exact) {
  avl_node_t *node;
  avl_compare_t cmp;
  void *userdata;
  int c;

  if (!exact)
    exact = &c;

  if (!tree) {
    *exact = 0;
    return nullptr;
  }

  node = tree->top;
  if (!node) {
    *exact = 0;
    return nullptr;
  }

  cmp = tree->cmpitem;
  userdata = tree->userdata;

  for (;;) {
    c = cmp(item, node->item, userdata);

    if (c < 0) {
      if (node->left)
        node = node->left;
      else
        return *exact = 0, node;
    } else if (c > 0) {
      if (node->right)
        node = node->right;
      else
        return *exact = 0, node->next;
    } else {
      return *exact = 1, node;
    }
  }

  return nullptr; /* To silence the compiler */
}
#endif

/* Searches for an item, returning either some exact
 * match, or (if no exact match could be found) the last (rightmost)
 * of the nodes that have an item smaller than the search item.
 * If exact is not nullptr, *exact will be set to:
 *    0  if the returned node is unequal or nullptr
 *    1  if the returned node is equal
 * Returns nullptr if no equal or smaller element could be found.
 * O(lg n) */
static AvlNode *AvlSearchRightish(const AvlTree *tree, const void *item, int *exact) {
  AvlNode *node;
  AvlCompare cmp;
  void *userdata;
  int c;

  if (!exact)
    exact = &c;

  if (!tree) {
    *exact = 0;
    return nullptr;
  }

  node = tree->top;
  if (!node) {
    *exact = 0;
    return nullptr;
  }

  cmp = tree->cmpitem;
  userdata = tree->userdata;

  for (;;) {
    c = cmp(item, node->item, userdata);

    if (c < 0) {
      if (node->left)
        node = node->left;
      else
        return *exact = 0, node->prev;
    } else if (c > 0) {
      if (node->right)
        node = node->right;
      else
        return *exact = 0, node;
    } else {
      return *exact = 1, node;
    }
  }

  return nullptr; /* To silence the compiler */
}

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_item_search_left(const avl_tree_t *tree, const void *item, int *exact) {
  avl_node_t *node;
  int c;

  if (!exact)
    exact = &c;

  if (!tree) {
    *exact = 0;
    return nullptr;
  }

  node = avl_search_leftish(tree, item, exact);
  if (*exact)
    return const_cast<avl_node_t *>(avl_search_leftmost_equal(tree, node, item));

  return const_cast<avl_node_t *>(node);
}
#endif

/* Searches for an item, returning either the last (rightmost) exact
 * match, or (if no exact match could be found) the last (rightmost)
 * of the nodes that have an item smaller than the search item.
 * If exact is not nullptr, *exact will be set to:
 *    0  if the returned node is inequal or nullptr
 *    1  if the returned node is equal
 * Returns nullptr if no equal or smaller element could be found.
 * O(lg n) */
static AvlNode *AvlItemSearchRight(const AvlTree *tree, const void *item, int *exact) {
  const AvlNode *node;
  int c;

  if (!exact)
    exact = &c;

  node = AvlSearchRightish(tree, item, exact);
  if (*exact)
    return const_cast<AvlNode *>(AvlSearchRightmostEqual(tree, node, item));

  return const_cast<AvlNode *>(node);
}

/* Searches for the item in the tree and returns a matching node if found
 * or nullptr if not.
 * O(lg n) */
static AvlNode *AvlItemSearch(const AvlTree *avltree, const void *item) {
  int c;
  AvlNode *n;
  n = AvlSearchRightish(avltree, item, &c);
  return c ? n : nullptr;
}

/* Initializes a new tree for elements that will be ordered using
 * the supplied strcmp()-like function.
 * Returns the value of avltree (even if it's nullptr).
 * O(1) */
static AvlTree *AvlTreeInit(AvlTree *avltree, AvlCompare cmp, AvlFree free_item) {
  if (avltree) {
    avltree->head = nullptr;
    avltree->tail = nullptr;
    avltree->top = nullptr;
    avltree->cmpitem = cmp;
    avltree->freeitem = free_item;
    avltree->userdata = nullptr;
  }
  return avltree;
}

/* Allocates and initializes a new tree for elements that will be
 * ordered using the supplied strcmp()-like function.
 * Returns nullptr if memory could not be allocated.
 * O(1) */
static AvlTree *AvlTreeConstruct(AvlCompare cmp, AvlFree free_item) { return AvlTreeInit(new AvlTree, cmp, free_item); }

/* Reinitializes the tree structure for reuse. Nothing is free()d.
 * Compare and free functions are left alone.
 * Returns the value of avltree (even if it's nullptr).
 * O(1) */
static AvlTree *AvlTreeClear(AvlTree *avltree) {
  if (avltree)
    avltree->top = avltree->head = avltree->tail = nullptr;
  return avltree;
}

static void AvlNodeFree(AvlTree *, AvlNode *node) { delete node; }

/* Free()s all nodes in the tree but leaves the tree itself.
 * If the tree's free is not nullptr it will be invoked on every item.
 * Returns the value of avltree (even if it's nullptr).
 * O(n) */
static AvlTree *AvlTreePurge(AvlTree *avltree) {
  AvlNode *node, *next;
  AvlFree func;
  void *userdata;

  if (!avltree)
    return nullptr;

  userdata = avltree->userdata;

  func = avltree->freeitem;

  for (node = avltree->head; node; node = next) {
    next = node->next;
    if (func)
      func(node->item, userdata);
    delete node;
  }

  return AvlTreeClear(avltree);
}

/* Frees the entire tree efficiently. Nodes will be free()d.
 * If the tree's free is not nullptr it will be invoked on every item.
 * O(n) */
static void AvlTreeDestruct(AvlTree *avltree) {
  if (!avltree)
    return;
  (void)AvlTreePurge(avltree);
  delete avltree;
}

/* Destroys an AVL tree when used as the deleter of a std::unique_ptr. */
struct AvlTreeDeleter {
  void operator()(AvlTree *avltree) const { AvlTreeDestruct(avltree); }
};

/* An AVL tree that is destroyed (including its items, if it has a free function) with its owner. */
using AvlTreePtr = std::unique_ptr<AvlTree, AvlTreeDeleter>;

static void AvlNodeClear(AvlNode *newnode) {
  newnode->left = newnode->right = nullptr;
#ifdef AVL_COUNT
  newnode->count = 1;
#endif
#ifdef AVL_DEPTH
  newnode->depth = 1;
#endif
}

/* Initializes memory for use as a node.
 * Returns the value of avlnode (even if it's nullptr).
 * O(1) */
static AvlNode *AvlNodeInit(AvlNode *newnode, const void *item) {
  if (newnode)
    newnode->item = const_cast<void *>(item);
  return newnode;
}

/* Allocates and initializes memory for use as a node.
 * Returns the value of avlnode (or nullptr if the allocation failed).
 * O(1) */
static AvlNode *AvlAlloc(AvlTree *, const void *item) { return AvlNodeInit(new AvlNode, item); }

/* Insert a node in an empty tree. If avl_node is nullptr, the tree will be
 * cleared and ready for re-use.
 * If the tree is not empty, the old nodes are left dangling.
 * O(1) */
static AvlNode *AvlInsertTop(AvlTree *avltree, AvlNode *newnode) {
  AvlNodeClear(newnode);
  newnode->prev = newnode->next = newnode->parent = nullptr;
  avltree->head = avltree->tail = avltree->top = newnode;
  return newnode;
}

/* Insert a node before another node. Returns the new node.
 * If old is nullptr, the item is appended to the tree.
 * O(lg n) */
static AvlNode *AvlNodeInsertBefore(AvlTree *avltree, AvlNode *node, AvlNode *newnode) {
  if (!avltree || !newnode)
    return nullptr;

  if (!node)
    return avltree->tail ? AvlNodeInsertAfter(avltree, avltree->tail, newnode) : AvlInsertTop(avltree, newnode);

  if (node->left)
    return AvlNodeInsertAfter(avltree, node->prev, newnode);

  AvlNodeClear(newnode);

  newnode->next = node;
  newnode->parent = node;

  newnode->prev = node->prev;
  if (node->prev)
    node->prev->next = newnode;
  else
    avltree->head = newnode;
  node->prev = newnode;

  node->left = newnode;
  AvlRebalance(avltree, node);
  return newnode;
}

/* Insert a node after another node. Returns the new node.
 * If old is nullptr, the item is prepended to the tree.
 * O(lg n) */
static AvlNode *AvlNodeInsertAfter(AvlTree *avltree, AvlNode *node, AvlNode *newnode) {
  if (!avltree || !newnode)
    return nullptr;

  if (!node)
    return avltree->head ? AvlNodeInsertBefore(avltree, avltree->head, newnode) : AvlInsertTop(avltree, newnode);

  if (node->right)
    return AvlNodeInsertBefore(avltree, node->next, newnode);

  AvlNodeClear(newnode);

  newnode->prev = node;
  newnode->parent = node;

  newnode->next = node->next;
  if (node->next)
    node->next->prev = newnode;
  else
    avltree->tail = newnode;
  node->next = newnode;

  node->right = newnode;
  AvlRebalance(avltree, node);
  return newnode;
}

/* Insert a node into the tree and return it.
 * Returns nullptr if an equal node is already in the tree.
 * O(lg n) */
static AvlNode *AvlNodeInsert(AvlTree *avltree, AvlNode *newnode) {
  AvlNode *node;
  int c;

  node = AvlSearchRightish(avltree, newnode->item, &c);
  return c ? nullptr : AvlNodeInsertAfter(avltree, node, newnode);
}

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_node_insert_left(avl_tree_t *avltree, avl_node_t *newnode) {
  return avl_node_insert_before(avltree, avl_item_search_left(avltree, newnode->item, nullptr), newnode);
}
#endif

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_node_insert_right(avl_tree_t *avltree, avl_node_t *newnode) {
  return avl_node_insert_after(avltree, avl_item_search_right(avltree, newnode->item, nullptr), newnode);
}
#endif

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_node_insert_somewhere(avl_tree_t *avltree, avl_node_t *newnode) {
  return avl_node_insert_after(avltree, avl_search_rightish(avltree, newnode->item, nullptr), newnode);
}
#endif

/* Insert an item into the tree and return the new node.
 * Returns nullptr and sets errno if memory for the new node could not be
 * allocated or if the node is already in the tree (EEXIST).
 * O(lg n) */
static AvlNode *AvlItemInsert(AvlTree *avltree, const void *item) {
  AvlNode *newnode;

  if (!avltree) {
    errno = EFAULT;
    return nullptr;
  }

  newnode = AvlAlloc(avltree, item);
  if (newnode) {
    if (AvlNodeInsert(avltree, newnode))
      return newnode;
    AvlNodeFree(avltree, newnode);
    errno = EEXIST;
  }
  return nullptr;
}

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_item_insert_somewhere(avl_tree_t *avltree, const void *item) {
  avl_node_t *newnode;

  if (!avltree) {
    errno = EFAULT;
    return nullptr;
  }

  newnode = avl_alloc(avltree, item);
  if (newnode)
    return avl_node_insert_somewhere(avltree, newnode);
  return nullptr;
}
#endif

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_item_insert_before(avl_tree_t *avltree, avl_node_t *node, const void *item) {
  avl_node_t *newnode;

  if (!avltree) {
    errno = EFAULT;
    return nullptr;
  }

  newnode = avl_alloc(avltree, item);
  if (newnode)
    return avl_node_insert_before(avltree, node, newnode);
  return nullptr;
}
#endif

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_item_insert_after(avl_tree_t *avltree, avl_node_t *node, const void *item) {
  avl_node_t *newnode;

  if (!avltree) {
    errno = EFAULT;
    return nullptr;
  }

  newnode = avl_alloc(avltree, item);
  if (newnode)
    return avl_node_insert_after(avltree, node, newnode);
  return nullptr;
}
#endif

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_item_insert_left(avl_tree_t *avltree, const void *item) {
  avl_node_t *newnode;

  if (!avltree) {
    errno = EFAULT;
    return nullptr;
  }

  newnode = avl_alloc(avltree, item);
  if (newnode)
    return avl_node_insert_left(avltree, newnode);
  return nullptr;
}
#endif

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_item_insert_right(avl_tree_t *avltree, const void *item) {
  avl_node_t *newnode;

  if (!avltree) {
    errno = EFAULT;
    return nullptr;
  }

  newnode = avl_alloc(avltree, item);
  if (newnode)
    return avl_node_insert_right(avltree, newnode);
  return nullptr;
}
#endif

/* Deletes a node from the tree.
 * Returns the value of the node (even if it's nullptr).
 * The item will NOT be free()d regardless of the tree's free handler.
 * This function comes in handy if you need to update the search key.
 * O(lg n) */
static AvlNode *AvlNodeUnlink(AvlTree *avltree, AvlNode *avlnode) {
  AvlNode *parent;
  AvlNode **superparent;
  AvlNode *subst, *left, *right;
  AvlNode *balnode;

  if (!avltree || !avlnode)
    return nullptr;

  if (avlnode->prev)
    avlnode->prev->next = avlnode->next;
  else
    avltree->head = avlnode->next;

  if (avlnode->next)
    avlnode->next->prev = avlnode->prev;
  else
    avltree->tail = avlnode->prev;

  parent = avlnode->parent;

  superparent = parent ? avlnode == parent->left ? &parent->left : &parent->right : &avltree->top;

  left = avlnode->left;
  right = avlnode->right;
  if (!left) {
    *superparent = right;
    if (right)
      right->parent = parent;
    balnode = parent;
  } else if (!right) {
    *superparent = left;
    left->parent = parent;
    balnode = parent;
  } else {
    subst = avlnode->prev;
    if (subst == left) {
      balnode = subst;
    } else {
      balnode = subst->parent;
      balnode->right = subst->left;
      if (balnode->right)
        balnode->right->parent = balnode;
      subst->left = left;
      left->parent = subst;
    }
    subst->right = right;
    subst->parent = parent;
    right->parent = subst;
    *superparent = subst;
  }

  AvlRebalance(avltree, balnode);

  return avlnode;
}

/* Deletes a node from the tree. Returns immediately if the node is nullptr.
 * If the tree's free is not nullptr, it is invoked on the item.
 * If it is, returns the item. In all other cases returns nullptr.
 * O(lg n) */
static void *AvlNodeDelete(AvlTree *avltree, AvlNode *avlnode) {
  void *item = nullptr;
  if (avlnode) {
    item = avlnode->item;
    (void)AvlNodeUnlink(avltree, avlnode);
    if (avltree->freeitem)
      avltree->freeitem(item, avltree->userdata);
    AvlNodeFree(avltree, avlnode);
  }
  return item;
}

/* Searches for an item in the tree and deletes it if found.
 * If the tree's free is not nullptr, it is invoked on the item.
 * If it is, returns the item. In all other cases returns nullptr.
 * O(lg n) */
static void *AvlItemDelete(AvlTree *avltree, const void *item) {
  return AvlNodeDelete(avltree, AvlItemSearch(avltree, item));
}

#if (!AVL_TREE_COMMENT_UNUSED)
static avl_node_t *avl_node_fixup(avl_tree_t *avltree, avl_node_t *newnode) {
  avl_node_t *oldnode = nullptr, *node;

  if (!avltree || !newnode)
    return nullptr;

  node = newnode->prev;
  if (node) {
    oldnode = node->next;
    node->next = newnode;
  } else {
    avltree->head = newnode;
  }

  node = newnode->next;
  if (node) {
    oldnode = node->prev;
    node->prev = newnode;
  } else {
    avltree->tail = newnode;
  }

  node = newnode->parent;
  if (node) {
    if (node->left == oldnode)
      node->left = newnode;
    else
      node->right = newnode;
  } else {
    oldnode = avltree->top;
    avltree->top = newnode;
  }

  return oldnode;
}
#endif

/**
 * AvlRebalance:
 * Rebalances the tree if one side becomes too heavy.  This function
 * assumes that both subtrees are AVL-trees with consistent data.  The
 * function has the additional side effect of recalculating the count of
 * the tree at this node.  It should be noted that at the return of this
 * function, if a rebalance takes place, the top of this subtree is no
 * longer going to be the same node.
 */
static void AvlRebalance(AvlTree *avltree, AvlNode *avlnode) {
  AvlNode *child;
  AvlNode *gchild;
  AvlNode *parent;
  AvlNode **superparent;

  while (avlnode) {
    parent = avlnode->parent;

    superparent = parent ? avlnode == parent->left ? &parent->left : &parent->right : &avltree->top;

    switch (AvlCheckBalance(avlnode)) {
    case -1:
      child = avlnode->left;
#ifdef AVL_DEPTH
      if (AVL_L_DEPTH(child) >= AVL_R_DEPTH(child)) {
#else
#ifdef AVL_COUNT
      if (AVL_L_COUNT(child) >= AVL_R_COUNT(child)) {
#else
#error No balancing possible.
#endif
#endif
        avlnode->left = child->right;
        if (avlnode->left)
          avlnode->left->parent = avlnode;
        child->right = avlnode;
        avlnode->parent = child;
        *superparent = child;
        child->parent = parent;
#ifdef AVL_COUNT
        avlnode->count = AVL_CALC_COUNT(avlnode);
        child->count = AVL_CALC_COUNT(child);
#endif
#ifdef AVL_DEPTH
        avlnode->depth = AVL_CALC_DEPTH(avlnode);
        child->depth = AVL_CALC_DEPTH(child);
#endif
      } else {
        gchild = child->right;
        avlnode->left = gchild->right;
        if (avlnode->left)
          avlnode->left->parent = avlnode;
        child->right = gchild->left;
        if (child->right)
          child->right->parent = child;
        gchild->right = avlnode;
        if (gchild->right)
          gchild->right->parent = gchild;
        gchild->left = child;
        if (gchild->left)
          gchild->left->parent = gchild;
        *superparent = gchild;
        gchild->parent = parent;
#ifdef AVL_COUNT
        avlnode->count = AVL_CALC_COUNT(avlnode);
        child->count = AVL_CALC_COUNT(child);
        gchild->count = AVL_CALC_COUNT(gchild);
#endif
#ifdef AVL_DEPTH
        avlnode->depth = AVL_CALC_DEPTH(avlnode);
        child->depth = AVL_CALC_DEPTH(child);
        gchild->depth = AVL_CALC_DEPTH(gchild);
#endif
      }
      break;
    case 1:
      child = avlnode->right;
#ifdef AVL_DEPTH
      if (AVL_R_DEPTH(child) >= AVL_L_DEPTH(child)) {
#else
#ifdef AVL_COUNT
      if (AVL_R_COUNT(child) >= AVL_L_COUNT(child)) {
#else
#error No balancing possible.
#endif
#endif
        avlnode->right = child->left;
        if (avlnode->right)
          avlnode->right->parent = avlnode;
        child->left = avlnode;
        avlnode->parent = child;
        *superparent = child;
        child->parent = parent;
#ifdef AVL_COUNT
        avlnode->count = AVL_CALC_COUNT(avlnode);
        child->count = AVL_CALC_COUNT(child);
#endif
#ifdef AVL_DEPTH
        avlnode->depth = AVL_CALC_DEPTH(avlnode);
        child->depth = AVL_CALC_DEPTH(child);
#endif
      } else {
        gchild = child->left;
        avlnode->right = gchild->left;
        if (avlnode->right)
          avlnode->right->parent = avlnode;
        child->left = gchild->right;
        if (child->left)
          child->left->parent = child;
        gchild->left = avlnode;
        if (gchild->left)
          gchild->left->parent = gchild;
        gchild->right = child;
        if (gchild->right)
          gchild->right->parent = gchild;
        *superparent = gchild;
        gchild->parent = parent;
#ifdef AVL_COUNT
        avlnode->count = AVL_CALC_COUNT(avlnode);
        child->count = AVL_CALC_COUNT(child);
        gchild->count = AVL_CALC_COUNT(gchild);
#endif
#ifdef AVL_DEPTH
        avlnode->depth = AVL_CALC_DEPTH(avlnode);
        child->depth = AVL_CALC_DEPTH(child);
        gchild->depth = AVL_CALC_DEPTH(gchild);
#endif
      }
      break;
    default:
#ifdef AVL_COUNT
      avlnode->count = AVL_CALC_COUNT(avlnode);
#endif
#ifdef AVL_DEPTH
      avlnode->depth = AVL_CALC_DEPTH(avlnode);
#endif
    }
    avlnode = parent;
  }
}

} // namespace coco
