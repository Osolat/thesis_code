//
// Created by benjamin on 2/11/22.
//

#include <cstdio>
#include <iostream>
#include <string>

#include "../lib/policy/policy_tree.h"

extern "C" {
#include <relic/relic.h>
}

#define NTESTS 5000

long long cpucycles(void) {
    unsigned long long result;
    asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
                 : "=a"(result)::"%rdx");
    return result;
}

static int cmp_llu(const void *a, const void *b) {
    if (*(unsigned long long *)a < *(unsigned long long *)b)
        return -1;
    if (*(unsigned long long *)a > *(unsigned long long *)b)
        return 1;
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
    uint32_t N_ATTR = atoi(argv[1]);

    core_init();

    bn_t order;
    bn_null(order);
    bn_new(order);

    pc_param_set_any();
    pc_param_print();
    pc_get_ord(order);

    bn_t rand;
    bn_null(rand);
    bn_new(rand);
    bn_rand_mod(rand, order);

    struct node tree_root;
    std::vector<policy_coefficient> res;
    tree_from_string(and_tree_formula(N_ATTR), &tree_root);
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();

        res = std::vector<policy_coefficient>();
        share_secret(&tree_root, rand, order, res, true);

    }
    printf("[");
    print_results("Results gen param():           ", t, NTESTS);

    bn_t attributes[N_ATTR];
    for (size_t i = 0; i < N_ATTR; i++) {
        bn_null(attributes[i]);
        bn_new(attributes[i]);
        bn_set_dig(attributes[i], i + 1);
    }

    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();

        try {
            check_satisfiability(&tree_root, attributes, N_ATTR);
        } catch (struct TreeUnsatisfiableException *e) {
            std::cout << e->what() << std::endl;
        }

        res = recover_coefficients(&tree_root, attributes, N_ATTR);
    }
    print_results("Results gen param():           ", t, NTESTS);
    printf("]\n");
    core_clean();
    return 0;
}
