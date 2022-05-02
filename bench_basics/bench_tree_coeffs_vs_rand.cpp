//
// Created by jonas on 4/1/22.
//

#include <cstdio>
#include <iostream>
#include <string>

#include "../lib/policy/policy_tree.h"

extern "C" {
#include <relic/relic.h>
}

#define NTESTS 3000

long long cpucycles(void) {
    unsigned long long result;
    asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
                 : "=a"(result)::"%rdx");
    return result;
}

static int cmp_llu(const void *a, const void *b) {
    if (*(unsigned long long *)a < *(unsigned long long *)b) return -1;
    if (*(unsigned long long *)a > *(unsigned long long *)b) return 1;
    return 0;
}

static unsigned long long median(unsigned long long *l, size_t llen) {
    qsort(l, llen, sizeof(unsigned long long), cmp_llu);

    if (llen % 2)
        return l[llen / 2];
    else
        return (l[llen / 2 - 1] + l[llen / 2]) / 2;
}

static unsigned long long average(unsigned long long *t, size_t tlen) {
    unsigned long long acc = 0;
    size_t i;
    for (i = 0; i < tlen; i++)
        acc += t[i];
    return acc / (tlen);
}

static void print_results(const char *s, unsigned long long *t, size_t tlen) {
    size_t i;
    for (i = 0; i < tlen - 1; i++) {
        t[i] = t[i + 1] - t[i];
    }
    printf("%llu,", average(t, tlen - 1));
}

unsigned long long t[NTESTS];

int main(int argc, char **argv) {
    if (argc == 1) {
        printf("Need to give argument\n");
        return 0;
    }

    int test_attr = atoi(argv[1]);
    uint32_t test_comp = test_attr;
    srand(time(NULL));

    core_init();
    bn_t order;
    pc_param_set_any();
    pc_param_print();
    g1_get_ord(order);

    g1_t g1_gen, g1_temp;
    g1_null(g1_gen);
    g1_new(g1_gen);
    g1_rand(g1_gen);
    g1_null(g1_temp);
    g1_new(g1_temp);

    gt_t gt_gen, gt_temp;
    gt_null(gt_gen);
    gt_new(gt_gen);
    gt_rand(gt_gen);
    gt_null(gt_temp);
    gt_new(gt_temp);

    g2_t g2_gen, g2_temp;
    g2_null(g2_gen);
    g2_new(g2_gen);
    g2_rand(g2_gen);
    g2_null(g2_temp);
    g2_new(g2_temp);

    bn_t rand_coeffs[test_attr];
    struct node tree_root;

    std::cout << "[g1 exp tree coeffs, g2 exp tree coeffs, gt exp tree coeffs, g1 exp rand coeffs, g2 exp rand coeffs, gt exp rand coeffs,]" << std::endl;

    for (int i = 0; i < test_comp; ++i) {
        bn_null(rand_coeffs[i]);
        bn_new(rand_coeffs[i]);
        bn_rand_mod(rand_coeffs[i], order);
    }

    tree_root = node();
    tree_from_string(and_tree_formula(test_attr), &tree_root);
    std::vector<policy_coefficient> res = std::vector<policy_coefficient>();
    share_secret(&tree_root, rand_coeffs[0], order, res, true);
    bn_t attributes[test_attr];

    for (size_t i = 0; i < test_attr; i++) {
        bn_null(attributes[i]);
        bn_new(attributes[i]);
        bn_set_dig(attributes[i], i + 1);
    }

    try {
        check_satisfiability(&tree_root, attributes, test_attr);
    } catch (struct TreeUnsatisfiableException *e) {
        std::cout << e->what() << std::endl;
    }

    res = recover_coefficients(&tree_root, attributes, test_attr);
    bn_t tree_coeffs[test_attr];
    for (auto it = res.begin(); it != res.end(); it++) {
        bn_null(tree_coeffs[it->leaf_index - 1]);
        bn_new(tree_coeffs[it->leaf_index - 1]);
        bn_copy(tree_coeffs[it->leaf_index - 1], it->coeff);
    }
    for (int i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        for (int j = 0; j < test_attr; ++j) {
            g1_mul(g1_temp, g1_gen, tree_coeffs[j]);
        }
    }
    printf("[");
    print_results("g1mulcoef:           ", t, NTESTS);

    for (int i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        for (int j = 0; j < test_attr; ++j) {
            g2_mul(g2_temp, g2_gen, tree_coeffs[j]);
        }
    }
    print_results("g2mulcoef:           ", t, NTESTS);
    
    for (int i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        for (int j = 0; j < test_attr; ++j) {
            gt_exp(gt_temp, gt_gen, tree_coeffs[j]);
        }
    }
    print_results("gtmulcoef:           ", t, NTESTS);
    

    for (int i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        for (int j = 0; j < test_attr; ++j) {
            g1_mul(g1_temp, g1_gen, rand_coeffs[j]);
        }
    }
    print_results("g1mulrand:           ", t, NTESTS);

    for (int i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        for (int j = 0; j < test_attr; ++j) {
            g2_mul(g2_temp, g2_gen, rand_coeffs[j]);
        }
    }
    print_results("g2mulrand:           ", t, NTESTS);
    
    for (int i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        for (int j = 0; j < test_attr; ++j) {
            gt_exp(gt_temp, gt_gen, rand_coeffs[j]);
        }
    }
    print_results("gtmulrand:           ", t, NTESTS);
    printf("]\n");
    core_clean();
    return 0;
}
