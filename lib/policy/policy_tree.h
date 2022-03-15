/*Created by Benjamin Birch Hansen on 2/22/22.*/

extern "C" {
#include <relic/relic.h>
}
#include <exception>
#include <string>
#include <vector>

enum gate_type { AND_GATE,
                 OR_GATE,
                 LEAF };

struct node {
    // TODO: Need to define what a node needs internally. Lagrange coefficients?
    enum gate_type gate;
    unsigned long long attribute_idx;
    bn_t attribute_zp;
    size_t children_num = 0;
    size_t leaf_index = 1;
    bn_t share;
    bn_t share_index;
    struct node *firstchild = NULL;
    struct node *nextsibling = NULL;
    bool marked_for_coeff = false;
};

struct policy_coefficient {
    // Index in a left-to-right indexing of leaves.
    size_t leaf_index;
    bn_t coeff;
    bn_t share;
};


/**
 * Initialises a tree from a boolean formula.
 * @param[in] formula			- "AND(OR(A), OR(B), OR(OR(C), OR(D))" = A AND B AND (C OR D)
 * @param[in] root              - Node to initialise as root.
 */
int tree_from_string(std::string formula, struct node *root);

/**
 * Add children to a node from formula. A result of parsing into a formula string.
 * @param[in] formula			- "OR(A), OR(B), OR(OR(C), OR(D)"
 * @param[in] parent            - parent node that we expect our additions to be childen of
 */
int add_children(struct node *parent, std::string formula);

/**
 * Prints entire tree
 * @param[in] root			- pointer to root of tree that we want printed
 */
void print_tree(struct node *root);

/**
 * Prints a single node in a tree
 * @param[in] n			- pointer to node we want printed
 */
int print_node(struct node *n);

std::string stringify_node(struct node *n);

/**
 * Share secret over access tree as in gpsw
 * @param[in] tree_root			- pointer to root of the access tree
 * @param[in] secret			- secret we want to share
 * @param[in] order			    - order of field
 * @param[in] res               - reference to vector to fill out with shares
 * @param[in] is_root           - indicates if this is the top of recursion
 */
int share_secret(struct node *tree_root, bn_t secret, bn_t order, std::vector<policy_coefficient> &res, bool is_root);

/**
 * Parses tree and checks if it can be satisfied by provided attributes. Also marks the minimal leafs
 * needed to find the coefficients that is needed to reconstruct the root secret
 * @param[in] tree_root			- pointer to root of the access tree
 * @param[in] attributes		- array of attributes
 * @param[in] num_attributes    - number of attributes
 * @throw                       - Throws exception if can't satisfy tree
 */
void check_satisfiability(struct node *tree_root, bn_t *attributes, size_t num_attributes);

/**
 * Parses tree and recovers needed coefficients to recover the secret at the root. Will first check if provided attributes can even satisfy policy.
 * @param[in] tree_root			- pointer to root of the access tree
 * @param[in] attributes		- array of attributes
 * @param[in] num_attributes    - number of attributes
 * @throw                       - Throws exception attributes do not satisfy tree
 */
std::vector<policy_coefficient> recover_coefficients(struct node *tree_root, bn_t *attributes, size_t num_attributes);

struct TreeUnsatisfiableException : public std::exception {
    const char *what() const throw() {
        return "Access policy tree could not be satisfied by attempted attributes";
    }
};

