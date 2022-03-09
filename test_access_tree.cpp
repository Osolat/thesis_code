#include <iostream>

//#include "lib/policy/policy_tree.h"
extern "C" {
#include "relic/relic.h"
}
/* int test1() {
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
} */

/* int test6() {
    bn_t secret, recovered, x[10], y[10];
    bn_null(secret);
    bn_new(secret);
    bn_null(recovered);
    bn_new(recovered);
    bn_t order;
    bn_null(order);
    bn_new(order);
    bn_set_dig(secret, 99);
    for (size_t i = 0; i < 10; i++) {
        bn_null(y[i]);
        bn_new(y[i]);
        bn_null(x[i]);
        bn_new(x[i]);
    }
    mpc_sss_gen(x, y, secret, order, 1, 10);
    mpc_sss_key(recovered, x, y, order, 1);
    bn_print(recovered);
} */

int main(int argc, char const *argv[]) {
    /*   std::cout << "Test 1 " << std::endl;
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
    bn_t secret, recovered, x[10], y[10];
    bn_null(secret);
    bn_new(secret);
    bn_null(recovered);
    bn_new(recovered);
    bn_t order;
    bn_null(order);
    bn_new(order);
    bn_set_dig(secret, 99);
    for (size_t i = 0; i < 10; i++) {
        bn_null(y[i]);
        bn_new(y[i]);
        bn_null(x[i]);
        bn_new(x[i]);
    }
    mpc_sss_gen(x, y, secret, order, 1, 10);
    mpc_sss_key(recovered, x, y, order, 1);
    //These are just to see if any function from relic_mpc compiles
    mpc_mt_gen(NULL, NULL);
    mpc_mt_lcl(NULL, NULL, NULL, NULL, NULL, NULL);
    mpc_mt_bct(NULL, NULL, NULL);
    mpc_mt_mul(NULL, NULL, NULL, NULL, NULL, NULL);
    md_map_sh512(NULL, NULL, NULL);
    rand_bytes(NULL,NULL);
    bn_print(recovered);
    return 0;
}