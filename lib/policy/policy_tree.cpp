#include "policy_tree.h"

#include <stdlib.h>

#include <algorithm>
#include <array>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

using namespace std;

const std::string WHITESPACE = " \n\r\t\f\v";
int id = 0;

std::string ltrim(const std::string& s) {
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}

std::string rtrim(const std::string& s) {
    size_t end = s.find_last_not_of(WHITESPACE);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

std::string trim(const std::string& s) {
    return rtrim(ltrim(s));
}

string find_first_gate(string formula) {
    string tokens[3] = {"AND", "OR", "LEAF"};
    int token_num = 3;
    size_t bestIndex = 0;
    size_t bestPosition = formula.size();
    string bestToken = "";
    for (size_t i = 0; i < token_num; ++i) {
        size_t const pos = formula.find(tokens[i]);
        if (pos >= bestPosition) {
            continue;
        }
        bestPosition = pos;
        bestIndex = i;
        bestToken = tokens[i];
    }
    return bestToken;
}

size_t find_index_first_gate(string formula) {
    string tokens[3] = {"AND", "OR", "LEAF"};
    int token_num = 3;
    size_t bestPosition = formula.size();
    for (size_t i = 0; i < token_num; ++i) {
        size_t const pos = formula.find(tokens[i]);
        if (pos >= bestPosition) {
            continue;
        }
        bestPosition = pos;
    }
    if (bestPosition == formula.size()) {
        bestPosition = -1;
    }
    return bestPosition;
}

size_t find_closing_paren(string text, size_t openPos) {
    size_t closePos = openPos;
    int counter = 1;
    while (counter > 0) {
        char c = text[++closePos];
        if (c == '(') {
            counter++;
        } else if (c == ')') {
            counter--;
        }
    }
    return closePos;
}

int tree_from_string(string formula, struct node* root) {
    // Absolutely needs formulas of form: "AND(OR(A), OR(B), OR(OR(C), OR(D))" = A AND B AND (C OR D)
    // So attributes are leaf OR nodes with 1 element
    // TODO: Do this in a better way than hardcoding node strings
    string bestToken = find_first_gate(formula);

    // TODO: push *root into function parameter
    // TODO: Initialise relic values
    if (bestToken.compare("AND") == 0) {
        root->gate = AND_GATE;
    } else {
        // TODO: Expand with more cases if we expand with other types of gates
        root->gate = OR_GATE;
    }
    // Strips AND(x) into x = arguments
    auto first = formula.find_first_of('(');
    auto last = formula.find_last_of(')');
    string arguments = formula.substr(first + 1, last - 1);
    if (arguments.substr(0, 5).find("attr") != string::npos) {
        // Formula was OR(attr...)
        // TODO: More stuff when we find a leaf node
        cout << "Tree with only a leaf? Was that intended?" << endl;
        root->gate = LEAF;
        return EXIT_SUCCESS;
    }

    add_children(root, arguments);

    return EXIT_SUCCESS;
}

int add_children(struct node* parent, string formula) {
    // Cases are either
    // 1: OR(attrx)
    // 2: AND/OR(something, something, something)
    size_t open = find_index_first_gate(formula);
    gate_type g;
    // Case 1
    if (formula.rfind("OR(attr", 0) == 0) {
        // TODO: More leaf additions
        g = LEAF;
        open = open + 2;
    } else {
        size_t found = formula.substr(open, open + 3).find("OR");
        // Case 2
        if (found != string::npos) {
            open = open + 2;
            g = OR_GATE;
        } else {
            found = formula.substr(open, open + 3).find("AND");
            if (found != string::npos) {
                open = open + 3;
                g = AND_GATE;
            } else {
                return 0;
            }
        }
    }
    cout << "Child add here" << endl;

    size_t closing = find_closing_paren(formula, open);
    struct node* leftmost_child = new node;
    leftmost_child->gate = g;
    leftmost_child->firstchild = nullptr;
    leftmost_child->nextsibling = nullptr;
    parent->firstchild = leftmost_child;
    if (g == OR_GATE || g == AND_GATE) {
        add_children(leftmost_child, formula.substr(open, closing - 1));
    }
    // Parse and add all siblings
    struct node* current_node = leftmost_child;
    string sibling_string = formula.substr(closing + 1, formula.size());
    sibling_string = trim(sibling_string);

    while (sibling_string.size() != 0 && sibling_string.at(0) == ',') {
        sibling_string = sibling_string.substr(1);
        sibling_string = trim(sibling_string);
        size_t open = find_index_first_gate(sibling_string);
        gate_type g;
        // Case 1
        if (sibling_string.rfind("OR(attr", 0) == 0) {
            // TODO: More leaf additions
            g = LEAF;
            open = open + 2;
        } else {
            size_t found = sibling_string.substr(open, open + 3).find("OR");
            // Case 2
            if (found != string::npos) {
                open = open + 2;
                g = OR_GATE;
            } else {
                found = sibling_string.substr(open, open + 3).find("AND");
                if (found != string::npos) {
                    open = open + 3;
                    g = AND_GATE;
                }
            }
        }
        cout << "Here sib" << endl;
        closing = find_closing_paren(sibling_string, open);
        struct node* brother = new node;
        brother->gate = g;
        brother->nextsibling = nullptr;
        brother->firstchild = nullptr;
        current_node->nextsibling = brother;
        if (g == OR_GATE || g == AND_GATE) {
            add_children(brother, sibling_string.substr(open + 1, closing - 3));
        }
        sibling_string = sibling_string.substr(closing + 1, sibling_string.size());
        current_node = brother;
    }
    return EXIT_SUCCESS;
}

int print_node(struct node* n) {
    string gate_name;
    switch (n->gate) {
        case LEAF:
            gate_name = "LEAF";
            break;
        case AND_GATE:
            gate_name = "AND";
            break;
        case OR_GATE:
            gate_name = "OR";
            break;
        default:
            cout << "Node has no gate_name!" << endl;
            return EXIT_FAILURE;
            break;
    }
    cout << gate_name << endl;
    return EXIT_SUCCESS;
}

void print_tree(struct node* root) {
    print_node(root);
    if (root->firstchild != nullptr) {
        cout << "child" << endl;
        print_tree(root->firstchild);
    }
    if (root->nextsibling != nullptr) {
        cout << "sibling" << endl;
        print_tree(root->nextsibling);
    }
}