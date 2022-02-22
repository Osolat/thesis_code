/*Created by Benjamin Birch Hansen on 2/22/22.*/
#include <relic/relic.h>

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
int tree_from_string(const char *formula);