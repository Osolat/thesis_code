#include "policy_tree.h"
#include <stdlib.h>

int tree_from_string(const char *formula) {
    //Absolutely needs formulas of form: "AND(OR(A), OR(B), OR(OR(C), OR(D))" = A AND B AND (C OR D)
    //So attributes are leaf OR nodes with 1 element
    struct node root;
    root.gate = AND_GATE;
    return EXIT_SUCCESS;
}
