/*Created by Benjamin Birch Hansen on 2/22/22.*/
#include <relic/relic.h>
#include <string>
enum gate_type{AND_GATE, OR_GATE};

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
 */
int tree_from_string(std::string formula);

/**
 * Add children to a node from formula. A result of parsing into a formula string.
 * @param[in] formula			- "OR(A), OR(B), OR(OR(C), OR(D)"
 * @param[in] parent            - parent node that we expect our additions to be childen of
 */
int add_children(struct node *parent, std::string formula);