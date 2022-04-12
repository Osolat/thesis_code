//
// Created by benjamin on 2/11/22.
//

#include <cstdio>
#include <iostream>
#include <string>

extern "C" {
#include <relic/relic.h>
}

#define NTESTS 10000

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

    g1_t g;
    g1_null(g);
    g1_new(g);
    g1_get_gen(g);

    g2_t h;
    g2_null(h);
    g2_new(h);
    g2_get_gen(h);

    g1_t pre_g[RLC_EP_TABLE_MAX];
    for (size_t i = 0; i < RLC_EP_TABLE_MAX; i++) {
        /* code */
        g1_null(pre_g[i]);
        g1_new(pre_g[i]);
    }
    g1_mul_pre(pre_g, g);

    g2_t pre_h[RLC_EP_TABLE_MAX];

    for (size_t i = 0; i < RLC_EP_TABLE_MAX; i++) {
        /* code */
        g2_null(pre_h[i]);
        g2_new(pre_h[i]);
    }
    g2_mul_pre(pre_h, h);

    bn_t bn_operands[NTESTS];

    for (size_t i = 0; i < NTESTS; i++) {
        bn_new(bn_operands[i]);
        bn_null(bn_operands[i]);
        bn_rand_mod(bn_operands[i], order);
    }

    g1_t temp;
    g1_new(temp);
    g1_null(temp);

    g2_t temp2;
    g2_new(temp2);
    g2_null(temp2);

    std::cout << "[mul_gen G1, mul_fix G1, mul_gen G2, mul_fix G2]" << std::endl;
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g1_mul_gen(temp, bn_operands[i]);
    }
    printf("[");
    print_results("Results gen param():           ", t, NTESTS);

    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g1_mul_fix(temp, pre_g, bn_operands[i]);
    }
    print_results("Results gen param():           ", t, NTESTS);

    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g2_mul_gen(temp2, bn_operands[i]);
    }
    
    print_results("Results gen param():           ", t, NTESTS);

    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g2_mul_fix(temp2, pre_h, bn_operands[i]);
    }
    print_results("Results gen param():           ", t, NTESTS);
    printf("]\n");

    core_clean();
    return 0;
}
