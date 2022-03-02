/*Created by Benjamin Birch Hansen on 2/22/22.*/
#include <relic/relic.h>
#include <string>
enum gate_type{AND_GATE, OR_GATE, LEAF};

struct node {
    // TODO: Need to define what a node needs internally. Lagrange coefficients?
    enum gate_type gate;
    bn_t share;
    struct node *firstchild;
    struct node *nextsibling;
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