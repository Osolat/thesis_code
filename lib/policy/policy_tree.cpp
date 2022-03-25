#include "policy_tree.h"

#include <stdlib.h>

#include <algorithm>
#include <array>
#include <iostream>
#include <regex>
#include <stack>
#include <string>
#include <vector>
using namespace std;

// TODO: This is very dangerous.
size_t global_leaf_idx;

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
    // string arguments = formula.substr(first + 1, last - 1);
    string arguments = string(&formula[first + 1], &formula[last]);
    if (arguments.substr(0, 5).find("attr") != string::npos) {
        // Formula was OR(attr...)
        // TODO: More stuff when we find a leaf node
        root->gate = LEAF;
        // Is LEAF gate
        string s = string(&formula[first + 5], &formula[last]);
        root->attribute_idx = stoull(s);
        bn_null(root->attribute_zp);
        bn_new(root->attribute_zp);
        bn_set_dig(root->attribute_zp, stoull(s));
        return EXIT_SUCCESS;
    }
    bn_null(root->share);
    bn_new(root->share);
    bn_zero(root->share);
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
    bn_null(leftmost_child->share);
    bn_new(leftmost_child->share);
    bn_zero(leftmost_child->share);
    if (g == OR_GATE || g == AND_GATE) {
        add_children(leftmost_child, string(&formula[open + 1], &formula[closing]));
    } else {
        // Is LEAF gate
        string s = string(&formula[open + 5], &formula[closing]);
        leftmost_child->attribute_idx = stoull(s);
        bn_null(leftmost_child->attribute_zp);
        bn_new(leftmost_child->attribute_zp);
        bn_set_dig(leftmost_child->attribute_zp, stoull(s));
    }
    // Parse and add all siblings
    struct node* current_node = leftmost_child;
    string sibling_string = string(&formula[closing + 1], &formula[formula.size()]);
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
        bn_null(brother->share);
        bn_new(brother->share);
        bn_zero(brother->share);
        if (g == OR_GATE || g == AND_GATE) {
            add_children(brother, string(&sibling_string[open + 1], &sibling_string[closing]));
        } else {
            // Is LEAF gate
            string s = string(&sibling_string[open + 5], &sibling_string[closing]);
            brother->attribute_idx = stoull(s);
            bn_null(brother->attribute_zp);
            bn_new(brother->attribute_zp);
            bn_set_dig(brother->attribute_zp, stoull(s));
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

int share_secret(struct node* tree_root, bn_t secret, bn_t order, std::vector<policy_coefficient>& res, bool root) {
    if (root) {
        global_leaf_idx = 1;
        if (tree_root->gate == LEAF) {
            policy_coefficient p = policy_coefficient();
            p.leaf_index = global_leaf_idx;
            tree_root->leaf_index = global_leaf_idx;
            bn_copy(tree_root->share, secret);
            bn_null(p.coeff);
            bn_new(p.coeff);
            bn_null(p.share);
            bn_new(p.share);
            bn_copy(p.share, secret);
            res.push_back(p);
            return 1;
        }
    }

    size_t children = tree_root->children_num;
    bn_t x[children], y[children];

    switch (tree_root->gate) {
        case AND_GATE: {
            //(n,n) treshold
            if (children >= 2) {
                /* code */

                for (size_t i = 0; i < children; i++) {
                    bn_null(x[i]);
                    bn_null(y[i]);
                    bn_new(x[i]);
                    bn_new(y[i]);
                }
                mpc_sss_gen(x, y, secret, order, children, children);
                struct node* child = tree_root->firstchild;

                int i = 0;
                while (child != NULL) {
                    bn_null(child->share);
                    bn_new(child->share);
                    bn_null(child->share_index);
                    bn_new(child->share_index);
                    bn_copy(child->share_index, x[i]);
                    bn_copy(child->share, y[i]);
                    share_secret(child, y[i], order, res, false);
                    child = child->nextsibling;
                    i++;
                }
                for (size_t i = 0; i < children; i++) {
                    bn_free(y[i]);
                    bn_free(x[i]);
                }
            } else {
                cout << "And node with 1 operand, fix your string." << endl;
            }
            break;
        }
        case OR_GATE: {
            //(1,n) threshold
            // Constant polynomial. f(0) = secret, so must be constant polynomial with coeff a_= secret.
            // All shares are the same
            struct node* child = tree_root->firstchild;
            bn_t bn_index;
            bn_null(bn_index);
            bn_new(bn_index);
            int i = 0;
            while (child != NULL) {
                bn_null(child->share);
                bn_new(child->share);
                bn_null(child->share_index);
                bn_new(child->share_index);
                bn_set_dig(bn_index, i);
                bn_copy(child->share_index, bn_index);
                bn_copy(child->share, secret);
                share_secret(child, secret, order, res, false);
                child = child->nextsibling;
                i++;
            }
            bn_free(bn_index);
            break;
        }
        case LEAF: {
            policy_coefficient p = policy_coefficient();
            p.leaf_index = global_leaf_idx;
            tree_root->leaf_index = global_leaf_idx;
            bn_null(p.coeff);
            bn_new(p.coeff);
            bn_null(p.share);
            bn_new(p.share);
            bn_copy(p.share, tree_root->share);
            res.push_back(p);
            global_leaf_idx++;
            break;
        }
        default:
            // Nothing needs to happen in a leaf.
            cout << "Reached undefined node in secret sharing function" << endl;
            break;
    }

    return 1;
}

int check_subtree_satisfiability(struct node* root, bn_t* attributes, size_t num_attributes) {
    switch (root->gate) {
        case AND_GATE: {
            /* code */
            // (n,n) threshold, all children must satisfy.
            struct node* child = root->firstchild;
            if (check_subtree_satisfiability(child, attributes, num_attributes) == 0) return 0;
            child = child->nextsibling;
            while (child != NULL) {
                if (check_subtree_satisfiability(child, attributes, num_attributes) == 0) return 0;
                child = child->nextsibling;
            }

            root->marked_for_coeff = true;
            return 1;
            break;
        }
        case LEAF: {
            for (size_t i = 0; i < num_attributes; i++) {
                if (bn_cmp(root->attribute_zp, attributes[i]) == 0) {
                    root->marked_for_coeff = true;
                    // Found correct attribute in leaf, so it can be satisfied
                    return 1;
                }
            }
            // No correct attribute in leaf, can't satisfy.
            return 0;
            break;
        }
        case OR_GATE: {
            // Here we must ensure only one branch is marked for efficiency.
            bool can_be_satisfied = false;
            struct node* child = root->firstchild;
            while (!can_be_satisfied && child != NULL) {
                can_be_satisfied = check_subtree_satisfiability(child, attributes, num_attributes);
                child = child->nextsibling;
            }
            if (can_be_satisfied) {
                root->marked_for_coeff = true;
                return 1;
            }
            return 0;
            break;
        }

        default:
            print_node(root);
            cout << "Found undefined gate_type. Everything will now break probably." << endl;
            return 0;
            break;
    }
}

void check_satisfiability(struct node* tree_root, bn_t* attributes, size_t num_attributes) {
    if (check_subtree_satisfiability(tree_root, attributes, num_attributes) == 0) throw TreeUnsatisfiableException();
}

std::vector<policy_coefficient> recover_coefficients(struct node* tree_root, bn_t* attributes, size_t num_attributes) {
    vector<policy_coefficient> result;
    std::stack<struct node*> node_stack;
    std::stack<bn_t*> coefficients;

    bn_t order;
    bn_null(order);
    bn_new(order);

    pc_get_ord(order);

    bn_t temp;
    bn_null(temp);
    bn_new(temp);

    if (num_attributes == 1) {
        policy_coefficient p = policy_coefficient();
        p.leaf_index = tree_root->leaf_index;
        bn_null(p.coeff);
        bn_new(p.coeff);
        bn_null(p.share);
        bn_new(p.share);

        bn_set_dig(temp, 1);
        bn_copy(p.coeff, temp);
        bn_copy(p.share, tree_root->share);
        result.push_back(p);
        return result;
    }

    bn_t unit;
    bn_null(unit);
    bn_new(unit);

    bn_t top_j, bot_j, bot_i, top_zero;
    bn_null(top_j);
    bn_new(top_j);
    bn_null(bot_j);
    bn_new(bot_j);
    bn_null(bot_i);
    bn_new(bot_i);
    bn_null(top_zero);
    bn_new(top_zero);
    bn_zero(top_zero);

    bn_set_dig(unit, 1);
    node_stack.push(tree_root);
    coefficients.push(&unit);

    struct node* current_node;
    struct node* child;
    struct node* brother;
    size_t threshold;
    while (!node_stack.empty()) {
        current_node = node_stack.top();
        bn_copy(temp, *coefficients.top());

        node_stack.pop();
        coefficients.pop();

        if (current_node->gate == LEAF) {
            policy_coefficient p = policy_coefficient();
            p.leaf_index = current_node->leaf_index;
            bn_null(p.coeff);
            bn_new(p.coeff);
            bn_null(p.share);
            bn_new(p.share);
            bn_copy(p.coeff, temp);
            bn_copy(p.share, current_node->share);
            result.push_back(p);
        } else {
            switch (current_node->gate) {
                case AND_GATE:
                    threshold = current_node->children_num;
                    break;
                case OR_GATE:
                    threshold = 1;
                    break;
                default:
                    cout << "Found Undefined node during coeff recovery, should not happen" << endl;
                    threshold = 0;
                    break;
            }
            size_t i = 1;
            bn_t* coeff = RLC_ALLOCA(bn_t, current_node->children_num);
            for (size_t i = 0; i < current_node->children_num; i++) {
                bn_null(coeff[i]);
                bn_new(coeff[i]);
                bn_set_dig(coeff[i], 1);
            }

            child = current_node->firstchild;
            while (child != NULL) {
                if (child->marked_for_coeff) {
                    size_t j = 1;
                    // Calculate langrange coefficients Prod(0 - j/ i - j)
                    brother = current_node->firstchild;
                    while (brother != NULL) {
                        if (brother->marked_for_coeff) {
                            if (j != i) {
                                // This seems to be a really slow way of doing Prod(0-j/i-j)
                                bn_set_dig(top_j, j);
                                bn_sub(top_j, top_zero, top_j);
                                bn_set_dig(bot_i, i);
                                bn_set_dig(bot_j, j);
                                bn_sub(bot_i, bot_i, bot_j);
                                bn_mod(bot_i, bot_i, order);
                                bn_mod_inv(bot_i, bot_i, order);

                                bn_mul(top_j, top_j, bot_i);

                                bn_mul(coeff[i - 1], coeff[i - 1], top_j);
                            }
                        }
                        j++;
                        brother = brother->nextsibling;
                    }
                    node_stack.push(child);
                    bn_mul(coeff[i - 1], temp, coeff[i - 1]);
                    coefficients.push(&coeff[i - 1]);
                }
                i++;
                child = child->nextsibling;
            }
        }
    }

    return result;
}

std::string and_tree_formula(size_t size) {
    string s;
    if (size == 1) {
        s.append("OR(attr1)");
        return s;
    }

    for (size_t i = 0; i < size - 2; i++) {
        string num = to_string(i + 1);
        s.append("AND(OR(attr");
        s.append(num);
        s.append("),");
    }

    string num = to_string(size - 1);
    s.append("AND(OR(attr");
    s.append(num);
    num = to_string(size);
    s.append("),OR(attr");
    s.append(num);

    for (size_t i = 0; i < size; i++) {
        s.append(")");
    }

    return s;
}

std::string or_tree_formula(size_t size) {
    string s;
    if (size == 1) {
        s.append("OR(attr1)");
        return s;
    }

    for (size_t i = 0; i < size - 2; i++) {
        string num = to_string(i + 1);
        s.append("OR(OR(attr");
        s.append(num);
        s.append("),");
    }

    string num = to_string(size - 1);
    s.append("OR(OR(attr");
    s.append(num);
    num = to_string(size);
    s.append("),OR(attr");
    s.append(num);

    for (size_t i = 0; i < size; i++) {
        s.append(")");
    }

    return s;
}