//
// Created by benjamin on 2/11/22.
//

#include <cstdio>
#include <iostream>
#include <string>

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
    int operands = atoi(argv[1]);

    core_init();

    bn_t order;
    bn_null(order);
    bn_new(order);

    pc_param_set_any();
    pc_param_print();
    pc_get_ord(order);

    gt_t m;
    gt_new(m);
    gt_null(m);

    g1_t g1_operands[NTESTS];
    g2_t g2_operands[NTESTS];
    gt_t gt_operands[NTESTS];

    for (size_t i = 0; i < NTESTS; i++) {
        g1_new(g1_operands[i]);
        g1_null(g1_operands[i]);
        g1_rand(g1_operands[i]);

        g2_new(g2_operands[i]);
        g2_null(g2_operands[i]);
        g2_rand(g2_operands[i]);

        gt_new(gt_operands[i]);
        gt_null(gt_operands[i]);
        gt_rand(gt_operands[i]);
    }

    g1_t temp;
    g1_new(temp);
    g1_null(temp);

    std::cout << "[negation g1, negate g2, negate(inv) gt]" << std::endl;
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g1_neg(g1_operands[i], g1_operands[i]);
    }
    printf("[");
    print_results("Results gen param():           ", t, NTESTS);

    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g2_neg(g2_operands[i], g2_operands[i]);
    }
    print_results("Results gen param():           ", t, NTESTS);

    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        gt_inv(gt_operands[i], gt_operands[i]);
    }
    print_results("Results gen param():           ", t, NTESTS);
    printf("]\n");

    core_clean();
    return 0;
}
