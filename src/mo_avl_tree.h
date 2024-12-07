#pragma once

/* We need either depths, counts or both (the latter being the default) */
#define AVL_DEPTH
#define AVL_COUNT

/* User supplied function to compare two items like strcmp() does.
 * For example: compare(a,b) will return:
 *   -1  if a < b
 *    0  if a = b
 *    1  if a > b
 */
typedef int (*avl_compare_t)(const void *a, const void *b, void *userdata);

/* User supplied function to delete an item when a node is free()d.
 * If NULL, the item is not free()d.
 */
typedef void (*avl_free_t)(void *item, void *userdata);

typedef struct avl_node {
  struct avl_node *prev;
  struct avl_node *next;
  struct avl_node *parent;
  struct avl_node *left;
  struct avl_node *right;
  void *item;
#ifdef AVL_COUNT
  unsigned long count;
#endif
#ifdef AVL_DEPTH
  unsigned char depth;
#endif
} avl_node_t;

typedef struct avl_tree {
  avl_node_t *top;
  avl_node_t *head;
  avl_node_t *tail;
  avl_compare_t cmpitem;
  avl_free_t freeitem;
  void *userdata;
  struct avl_allocator *allocator;
} avl_tree_t;

avl_node_t *avl_alloc(avl_tree_t *avltree,const void *item);
avl_node_t *avl_insert_top(avl_tree_t *avltree,avl_node_t *newnode);
avl_node_t *avl_item_insert(avl_tree_t *avltree,const void *item);
avl_node_t *avl_item_search(const avl_tree_t *avltree,const void *item);
avl_node_t *avl_item_search_right(const avl_tree_t *tree,const void *item,int *exact);
avl_node_t *avl_node_init(avl_node_t *newnode,const void *item);
avl_node_t *avl_node_insert(avl_tree_t *avltree,avl_node_t *newnode);
avl_node_t *avl_node_insert_after(avl_tree_t *avltree,avl_node_t *node,avl_node_t *newnode);
avl_node_t *avl_node_insert_before(avl_tree_t *avltree,avl_node_t *node,avl_node_t *newnode);
avl_node_t *avl_node_unlink(avl_tree_t *avltree,avl_node_t *avlnode);
avl_node_t *avl_search_rightish(const avl_tree_t *tree,const void *item,int *exact);
avl_tree_t *avl_tree_clear(avl_tree_t *avltree);
avl_tree_t *avl_tree_construct(avl_compare_t cmp,avl_free_t free_item);
avl_tree_t *avl_tree_init(avl_tree_t *avltree,avl_compare_t cmp,avl_free_t free_item);
avl_tree_t *avl_tree_purge(avl_tree_t *avltree);
const avl_node_t *avl_search_rightmost_equal(const avl_tree_t *tree,const avl_node_t *node,const void *item);
int avl_check_balance(avl_node_t *avlnode);
void *avl_item_delete(avl_tree_t *avltree,const void *item);
void *avl_node_delete(avl_tree_t *avltree,avl_node_t *avlnode);
void avl_node_clear(avl_node_t *newnode);
void avl_node_free(avl_tree_t *avltree,avl_node_t *node);
void avl_rebalance(avl_tree_t *,avl_node_t *);
void avl_tree_destruct(avl_tree_t *avltree);
unsigned long avl_count(const avl_tree_t *avltree);
