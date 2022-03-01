#include "lib/policy/policy_tree.h"

int test1() {
    char formula[] = "AND(OR(attr1), OR(attr2), OR(OR(attr3), OR(attr4))";
    struct node root;
    tree_from_string(formula, &root);
    return EXIT_SUCCESS;
}

int main(int argc, char const *argv[]) {
    test1();
    return 0;
}