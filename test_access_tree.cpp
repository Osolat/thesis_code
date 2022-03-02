#include <iostream>

#include "lib/policy/policy_tree.h"

int test1() {
    char formula[] = "AND(OR(attr1),OR(attr2),OR(OR(attr3),OR(attr4))";
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
    char formula[] = "AND(OR(attr1),OR(attr2),OR(OR(attr3),OR(attr4))";
    struct node root;
    tree_from_string(formula, &root);
    std::cout << formula << std::endl;
    print_tree(&root);
    return EXIT_SUCCESS;
}

int main(int argc, char const *argv[]) {
    /*   struct node root;
      root.gate = AND_GATE;
      struct node* current_node = new node;
      root.firstchild = current_node;
      current_node -> gate = LEAF;

      for (size_t i = 0; i < 5; i++)
      {

          struct node* brother = new node;
          brother->gate = LEAF;
          current_node->nextsibling = brother;
          current_node = brother;
      }
      print_tree(&root); */
    std::cout << "Test 1 " << std::endl;
    test1();
    std::cout << "--------------------- " << std::endl;
    std::cout << "Test 2 " << std::endl;
    test2();
    std::cout << "--------------------- " << std::endl;
    std::cout << "Test 3 " << std::endl;
    test3();
    return 0;
}