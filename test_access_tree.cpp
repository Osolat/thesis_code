#include <iostream>

#include "lib/policy/policy_tree.h"

extern "C" {
#include <relic/relic_test.h>
}
int test1() {
    char formula[] = "AND(OR(attr1),OR(attr2),OR(OR(attr3),OR(attr4)))";
    struct node root = node();
    tree_from_string(formula, &root);
    return EXIT_SUCCESS;
}

int test2() {
    struct node n = node();
    n.gate = LEAF;
    struct node n2 = node();
    n2.gate = AND_GATE;
    struct node n3 = node();
    n3.gate = OR_GATE;
    print_node(&n);
    print_node(&n2);
    print_node(&n3);
    return EXIT_SUCCESS;
}

int test3() {
    char formula[] = "AND(OR(attr1),OR(attr2),OR(OR(attr3),OR(attr4)))";
    struct node root;
    tree_from_string(formula, &root);
    std::cout << formula << std::endl;
    print_tree(&root);
    return EXIT_SUCCESS;
}

int test4() {
    char formula[] = "AND(OR(attr1),OR(attr2),OR(OR(attr3),AND(OR(attr4),OR(attr5))))";
    struct node root;
    tree_from_string(formula, &root);
    std::cout << formula << std::endl;
    print_tree(&root);
    return EXIT_SUCCESS;
}

int test5() {
    char formula[] = "AND(OR(attr1),OR(attr2),OR(attr3),OR(attr4),OR(attr5),OR(attr6),OR(attr7),OR(attr8),OR(attr9),OR(attr10),OR(attr11))";
    struct node root;
    tree_from_string(formula, &root);
    std::cout << formula << std::endl;
    print_tree(&root);
    return EXIT_SUCCESS;
}

int test6() {
    core_init();
    char formula[] = "AND(OR(attr1),OR(attr2),OR(OR(attr3),OR(attr4)))";
    struct node root;
    tree_from_string(formula, &root);

    bn_t secret, order;
    bn_null(secret);
    bn_new(secret);
    bn_null(order);
    bn_new(order);
    pc_param_set_any();
    g1_get_ord(order);
    bn_set_dig(secret, 99);
    share_secret(&root, secret, order);
    print_tree(&root);

    bn_t attributes[4];
    for (size_t i = 0; i < 4; i++) {
        bn_null(attributes[i]);
        bn_new(attributes[i]);
        bn_set_dig(attributes[i], i + 1);
    }
    bn_t fail_attributes[3];
    for (size_t i = 0; i < 3; i++) {
        bn_null(fail_attributes[i]);
        bn_new(fail_attributes[i]);
        bn_set_dig(fail_attributes[i], i + 2);
    }

    try {
        check_satisfiability(&root, attributes, 4);
        std::cout << "Satisfiable with correct attributes" << std::endl;
    } catch (struct TreeUnsatisfiableException *e) {
        std::cout << e->what() << std::endl;
    }

    try {
        check_satisfiability(&root, fail_attributes, 3);
    } catch (struct TreeUnsatisfiableException *e) {
        std::cout << e->what() << std::endl;
    }
    core_clean();
}

int main(int argc, char const *argv[]) {
    core_init();

    /* std::cout << "Test 1 " << std::endl;
    test1();
    std::cout << "--------------------- " << std::endl;
    std::cout << "Test 2 " << std::endl;
    test2();
    std::cout << "--------------------- " << std::endl;
    std::cout << "Test 3 " << std::endl;
    test3();
    std::cout << "--------------------- " << std::endl;
    std::cout << "Test 4 " << std::endl;
    test4();
    std::cout << "--------------------- " << std::endl;
    std::cout << "Test 5 " << std::endl;
    test5(); */
    core_init();
    char formula[] = "AND(OR(attr1),OR(attr2),OR(attr3))";
    struct node root;
    tree_from_string(formula, &root);

    bn_t secret, order;
    bn_null(secret);
    bn_new(secret);
    bn_null(order);
    bn_new(order);
    pc_param_set_any();
    g1_get_ord(order);
    bn_set_dig(secret, 99);
    bn_print(secret);
    share_secret(&root, secret, order);
    print_tree(&root);

    bn_t attributes[3];
    for (size_t i = 0; i < 3; i++) {
        bn_null(attributes[i]);
        bn_new(attributes[i]);
        bn_set_dig(attributes[i], i + 1);
    }

    try {
        check_satisfiability(&root, attributes, 3);
        std::cout << "Satisfiable with correct attributes" << std::endl;
    } catch (struct TreeUnsatisfiableException *e) {
        std::cout << e->what() << std::endl;
    }
    std::vector<policy_coefficient> res;
    res = recover_coefficients(&root, attributes, 3);
    bn_t gather;
    bn_null(gather);
    bn_new(gather);
    bn_zero(gather);

    bn_t tmp;
    bn_null(tmp);
    bn_new(tmp);

    for (std::vector<policy_coefficient>::iterator it = res.begin(); it != res.end(); ++it) {
        bn_print(it->share);
        bn_print(it->coeff);
        bn_mul(tmp, it->coeff, it->share);
        bn_add(gather, gather, tmp);
        bn_mod(gather, gather, order);
    }
    bn_print(secret);
    bn_print(gather);
    core_clean();

    return 1;
}