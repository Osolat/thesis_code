#include <iostream>

#include "lib/policy/policy_tree.h"

extern "C" {
#include <relic/relic_test.h>
}

#include "legacy/arith/l_zelement_bp.h"
#include "legacy/arith/l_zgroup.h"
#include "legacy/l_zattributelist.h"
#include "legacy/l_zbytestring.h"
#include "legacy/l_zconstants.h"
#include "legacy/l_zfunctioninput.h"
#include "legacy/l_zinteger.h"
#include "legacy/l_zobject.h"
#include "legacy/l_zpolicy.h"
#include "legacy/lsss/l_zlsss.h"

extern "C" {
#include "legacy/arith/l_zelement.h"
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
    // char formula[] = "AND(AND(OR(attr1),OR(attr2)),AND(OR(attr3),OR(attr4),OR(attr5)))";
    //char formula[] = "AND(OR(attr1),OR(attr2),OR(OR(attr3),OR(attr4),OR(attr5)))";

    bn_t two;
    bn_null(two);
    bn_new(two);

    bn_t nottwo;
    bn_null(nottwo);
    bn_new(nottwo);

    bn_set_dig(two, 2);
    bn_neg(nottwo, two);
    bn_print(two);
    bn_print(nottwo);
    char formula[] = "AND(OR(attr1),OR(attr2),OR(attr3),OR(attr4),OR(attr5))";
    struct node root;
    tree_from_string(formula, &root);

    bn_t secret, order;
    bn_null(secret);
    bn_new(secret);
    bn_null(order);
    bn_new(order);
    pc_param_set_any();
    g1_get_ord(order);
    bn_set_dig(secret, 105);
    bn_print(secret);
    std::vector<policy_coefficient> res;
    share_secret(&root, secret, order, res, true);
    print_tree(&root);

    for (std::vector<policy_coefficient>::iterator it = res.begin(); it != res.end(); ++it) {
        bn_print(it->share);
        std::cout << it->leaf_index << std::endl;
    }

    bn_t attributes[5];
    for (size_t i = 0; i < 5; i++) {
        bn_null(attributes[i]);
        bn_new(attributes[i]);
        bn_set_dig(attributes[i], i + 1);
    }

    try {
        check_satisfiability(&root, attributes, 5);
        std::cout << "Satisfiable with correct attributes" << std::endl;
    } catch (struct TreeUnsatisfiableException *e) {
        std::cout << e->what() << std::endl;
    }
    res = std::vector<policy_coefficient>();
    res = recover_coefficients(&root, attributes, 5);
    bn_t gather;
    bn_null(gather);
    bn_new(gather);
    bn_zero(gather);

    bn_t tmp;
    bn_null(tmp);
    bn_new(tmp);

    for (std::vector<policy_coefficient>::iterator it = res.begin(); it != res.end(); ++it) {
        //bn_print(it->share);
        bn_print(it->coeff);
        //std::cout << it->leaf_index << std::endl;
        bn_mul(tmp, it->coeff, it->share);
        bn_add(gather, gather, tmp);
    }
    bn_mod(gather, gather, order);
    bn_print(secret);
    bn_print(gather);
    return EXIT_SUCCESS;
}

int main(int argc, char const *argv[]) {
    core_init();

    /*  std::cout << "Test 1 " << std::endl;
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
     test5();
     std::cout << "--------------------- " << std::endl;
     std::cout << "Test 6 " << std::endl;
     test6(); */

    test6();
    std::cout << and_tree_formula(10) << std::endl;
    struct node root;
    std::string str = "AND(AND(AND(AND(OR(attr5), OR(attr4)), OR(attr3)), OR(attr2)), OR(attr1))";
    tree_from_string(str, &root);
    bn_t secret, order;
    bn_null(secret);
    bn_new(secret);
    bn_null(order);
    bn_new(order);
    pc_param_set_any();
    g1_get_ord(order);
    bn_print(order);
    bn_set_dig(secret, 105);
    bn_print(secret);
    std::vector<policy_coefficient> res;
    share_secret(&root, secret, order, res, true);
    // print_tree(&root);

    for (std::vector<policy_coefficient>::iterator it = res.begin(); it != res.end(); ++it) {
        bn_print(it->share);
        std::cout << it->leaf_index << std::endl;
    }

    bn_t attributes[5];
    for (size_t i = 0; i < 5; i++) {
        bn_null(attributes[i]);
        bn_new(attributes[i]);
        bn_set_dig(attributes[i], i + 1);
    }

    try {
        check_satisfiability(&root, attributes, 5);
        std::cout << "Satisfiable with correct attributes" << std::endl;
    } catch (struct TreeUnsatisfiableException *e) {
        std::cout << e->what() << std::endl;
    }
    res = std::vector<policy_coefficient>();
    res = recover_coefficients(&root, attributes, 5);
    bn_t gather;
    bn_null(gather);
    bn_new(gather);
    bn_zero(gather);

    bn_t tmp;
    bn_null(tmp);
    bn_new(tmp);
    std::cout << "Boy Im here" << std::endl;
    for (std::vector<policy_coefficient>::iterator it = res.begin(); it != res.end(); ++it) {
        bn_print(it->share);
        bn_print(it->coeff);
        std::cout << it->leaf_index << std::endl;
        bn_mul(tmp, it->coeff, it->share);
        bn_add(gather, gather, tmp);
    }
    bn_mod(gather, gather, order);
    bn_print(secret);
    bn_print(gather);
    core_clean();

    return 1;
}