/*Created by Benjamin Birch Hansen on 2/22/22.*/

extern "C" {
#include <relic/relic.h>
}
#include <string>
enum gate_type { AND_GATE,
                 OR_GATE,
                 LEAF };

struct node {
    // TODO: Need to define what a node needs internally. Lagrange coefficients?
    enum gate_type gate;
    unsigned long long attribute_idx;
    size_t children_num = 0;
    bn_t share;
    bn_t share_index;
    struct node *firstchild = NULL;
    struct node *nextsibling = NULL;
};

/**
 * Initialises a tree from a boolean formula.
 * @param[in] formula			- "AND(OR(A), OR(B), OR(OR(C), OR(D))" = A AND B AND (C OR D)
 * @param[in] root              - Node to initialise as root.
 */
int tree_from_string(std::string formula, struct node *root);

/**
 * Add children to a node from formula. A result of parsing into a formula string.
 * @param[in] formula			- "OR(A), OR(B), OR(OR(C), OR(D)"
 * @param[in] parent            - parent node that we expect our additions to be childen of
 */
int add_children(struct node *parent, std::string formula);

/**
 * Prints entire tree
 * @param[in] root			- pointer to root of tree that we want printed
 */
void print_tree(struct node *root);

/**
 * Prints a single node in a tree
 * @param[in] n			- pointer to node we want printed
 */
int print_node(struct node *n);

std::string stringify_node(struct node *n);

/**
 * Share secret over access tree as in gpsw
 * @param[in] tree_root			- pointer to root of the access tree
 * @param[in] secret			- secret we want to share
 * @param[in] order			    - order of field
 */
int share_secret(struct node *tree_root, bn_t secret, bn_t order);