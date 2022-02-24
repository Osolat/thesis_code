#include "policy_tree.h"

#include <stdlib.h>

#include <array>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

int tree_from_string(string formula) {
    // Absolutely needs formulas of form: "AND(OR(A), OR(B), OR(OR(C), OR(D))" = A AND B AND (C OR D)
    // So attributes are leaf OR nodes with 1 element
    // TODO: Do this in a better way than hardcoding node strings

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
    // TODO push *root into function parameter
    struct node root;
    if (bestToken.compare("AND") == 0) {
        root.gate = AND_GATE;
    } else {
        // TODO: Expand with more cases if we expand with other types of gates
        root.gate = OR_GATE;
    }
    // Strips AND(x) into x = arguments
    auto first = formula.find_first_of('(');
    auto last = formula.find_last_of(')');
    string arguments = formula.substr(first + 1, last + 1);

    add_children(&root, arguments);

    return EXIT_SUCCESS;
}

int add_children(struct node* parent, string formula) {
    stringstream ss(formula);

    vector<string> arguments;
    while (ss.good()) {
        string substr;
        getline(ss, substr, ',');
        if (substr.find_first_not_of(' ') != std::string::npos) {
            // Not just a whitespace token
            printf("%s \n", substr.c_str());
            arguments.push_back(substr);
        }
    }
    return EXIT_SUCCESS;
}
