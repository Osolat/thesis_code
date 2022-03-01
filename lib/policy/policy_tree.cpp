#include "policy_tree.h"

#include <stdlib.h>

#include <array>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

using namespace std;

string find_first_gate(string formula) {
    string tokens[2] = {"AND", "OR"};
    int token_num = 2;
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
    string tokens[2] = {"AND", "OR"};
    int token_num = 2;
    size_t bestPosition = formula.size();
    for (size_t i = 0; i < token_num; ++i) {
        size_t const pos = formula.find(tokens[i]);
        if (pos >= bestPosition) {
            continue;
        }
        bestPosition = pos;
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
    return closePos - 1;
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
    string arguments = formula.substr(first + 1, last + 1);

    add_children(root, arguments);

    return EXIT_SUCCESS;
}

int add_children(struct node* parent, string formula) {
    /*  stringstream ss(formula);

     vector<string> arguments;
     while (ss.good()) {
         string substr;
         getline(ss, substr, ',');
         if (substr.find_first_not_of(' ') != std::string::npos) {
             // Not just a whitespace token
             printf("%s \n", substr.c_str());
             arguments.push_back(substr);
         }
     } */
    printf("%s \n", formula.c_str());
    vector<string> arguments;
    regex r("[a-zA-Z0-9\\s]+\\([a-zA-Z0-9_@&:.# -]+\\)");
    smatch m;
    regex_search(formula, m, r);
    for (auto x : m) {
        // cout << x << "\n ";
    }

    // Cases are either
    // 1: OR(attrx)
    // 2: AND/OR(something, something, something)
    size_t open = find_index_first_gate(formula);
    size_t found = formula.substr(open, open + 3).find("OR");
    gate_type g;
    if (found != string::npos) {
        open = open + 2;
        g = OR_GATE;
    } else {
        found = formula.substr(open, open + 3).find("AND");
        if (found != string::npos) {
            open = open + 3;
            g = AND_GATE;
        }
    }
    size_t closing = find_closing_paren(formula, open);

    
    return EXIT_SUCCESS;
}
