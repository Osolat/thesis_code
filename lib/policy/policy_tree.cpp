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
        size_t found = string(&formula[open], &formula[open + 4]).find("OR");
        // Case 2
        if (found != string::npos) {
            open = open + 2;
            g = OR_GATE;
        } else {
            found = string(&formula[open], &formula[open + 4]).find("AND");
            if (found != string::npos) {
                open = open + 3;
                g = AND_GATE;
            } else {
                return 0;
            }
        }
    }

    size_t closing = find_closing_paren(formula, open);
    struct node* leftmost_child = new node;
    leftmost_child->gate = g;
    parent->firstchild = leftmost_child;
    parent->children_num++;
    if (g == OR_GATE || g == AND_GATE) {
        add_children(leftmost_child, string(&formula[open], &formula[closing]));
    } else {
        // Is LEAF gate
        string s = string(&formula[open + 5], &formula[closing]);
        leftmost_child->attribute_idx = stoull(s);
    }
    // Parse and add all siblings
    struct node* current_node = leftmost_child;
    string sibling_string = formula.substr(closing + 1, formula.size());

    while (sibling_string.size() != 0 && sibling_string.at(0) == ',') {
        sibling_string = sibling_string.substr(1);
        size_t open = find_index_first_gate(sibling_string);
        gate_type g;
        // Case 1
        if (sibling_string.rfind("OR(attr", 0) == 0) {
            // TODO: More leaf additions
            g = LEAF;
            open = open + 2;
        } else {
            size_t found = string(&sibling_string[open], &sibling_string[open + 4]).find("OR");
            // Case 2
            if (found != string::npos) {
                open = open + 2;
                g = OR_GATE;
            } else {
                found = string(&sibling_string[open], &sibling_string[open + 4]).find("AND");
                if (found != string::npos) {
                    open = open + 3;
                    g = AND_GATE;
                }
            }
        }
        closing = find_closing_paren(sibling_string, open);
        struct node* brother = new node;
        brother->gate = g;
        current_node->nextsibling = brother;
        parent->children_num++;
        if (g == OR_GATE || g == AND_GATE) {
            add_children(brother, string(&sibling_string[open + 1], &sibling_string[closing]));
        } else {
            // Is LEAF gate
            string s = string(&sibling_string[open + 5], &sibling_string[closing]);
            brother->attribute_idx = stoull(s);
        }
        sibling_string = sibling_string.substr(closing + 1);
        current_node = brother;
    }
    return EXIT_SUCCESS;
}

string stringify_node(struct node* n) {
    string s;
    switch (n->gate) {
        case LEAF:
            s = "LEAF";
            s.append(" - attribute").append(to_string(n->attribute_idx));
            break;
        case AND_GATE:
            s = "AND";
            break;
        case OR_GATE:
            s = "OR";
            break;
        default:
            s = "UNDEFINED NODE";
            break;
    }
    s.append(" with ").append(to_string(n->children_num)).append(" children");
    return s;
}

int print_node(struct node* n) {
    cout << stringify_node(n) << endl;
    return EXIT_SUCCESS;
}

void print_tree(struct node* root) {
    print_node(root);
    if (root->firstchild != NULL) {
        cout << "child" << endl;
        bn_print(root->firstchild->share);

        print_tree(root->firstchild);
    }
    if (root->nextsibling != NULL) {
        cout << "sibling" << endl;
        bn_print(root->nextsibling->share);

        print_tree(root->nextsibling);
    }
}

int share_secret(struct node* tree_root, bn_t secret, bn_t order) {
    size_t children = tree_root->children_num;
    bn_t x[children], y[children];

    // TODO This is because of Relic's Shamir implementation only supporting t >= 3
    switch (tree_root->gate) {
        case AND_GATE: {
            //(n,n) treshold
            if (children >= 3) {
                /* code */

                for (size_t i = 0; i < children; i++) {
                    bn_null(x[i]);
                    bn_null(y[i]);
                    bn_new(x[i]);
                    bn_new(y[i]);
                }
                mpc_sss_gen(x, y, secret, order, children, children);
                struct node* child = tree_root->firstchild;
                if (child != NULL) {
                    cout << "First child got shares" << endl;
                    bn_null(child->share);
                    bn_new(child->share);
                    bn_null(child->share_index);
                    bn_new(child->share_index);
                    bn_copy(child->share_index, x[0]);
                    bn_copy(child->share, y[0]);
                    share_secret(child, y[0], order);
                    // Get next child which is the first child's brother.
                    child = child->nextsibling;
                }
                int i = 1;
                while (child != NULL) {
                    bn_null(child->share);
                    bn_new(child->share);
                    bn_null(child->share_index);
                    bn_new(child->share_index);
                    bn_copy(child->share_index, x[i]);
                    bn_copy(child->share, y[i]);
                    share_secret(child, y[i], order);
                    child = child->nextsibling;
                    i++;
                }
                for (size_t i = 0; i < children; i++) {
                    bn_free(y[i]);
                    bn_free(x[i]);
                }
            }
            break;
        }
        case OR_GATE: {
            //(1,n) threshold
            // Constant polynomial. f(0) = secret, so must be constant polynomial with coeff a_= secret.
            // All shares are the same
            cout << "Or or or " << endl;
            struct node* child = tree_root->firstchild;
            bn_t bn_index;
            bn_null(bn_index);
            bn_new(bn_index);
            bn_set_dig(bn_index, 0);
            if (child != NULL) {
                cout << "Or not null " << endl;
                bn_null(child->share);
                bn_new(child->share);
                bn_null(child->share_index);
                bn_new(child->share_index);
                bn_copy(child->share_index, bn_index);
                bn_copy(child->share, secret);
                cout << "Or not null blyat " << endl;

                share_secret(child, secret, order);
                // Get next child which is the first child's brother.
                child = child->nextsibling;
            }
            int i = 1;
            while (child != NULL) {
                bn_null(child->share);
                bn_new(child->share);
                bn_null(child->share_index);
                bn_new(child->share_index);
                bn_set_dig(bn_index, i);
                bn_copy(child->share_index, bn_index);
                bn_copy(child->share, secret);
                share_secret(child, secret, order);
                child = child->nextsibling;
                i++;
            }
            bn_free(bn_index);
            break;
        }
        default:
            // Nothing needs to happen in a leaf.
            break;
    }

    return 1;
}