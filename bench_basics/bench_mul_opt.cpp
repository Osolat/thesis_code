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

    g1_t g1_operands[operands];
    g2_t g2_operands[operands];
    bn_t x[operands];
    std::cout << "[G1 iterative mul+add, G1 mul_sim_lot, G2 iterative mul+add, G2 mul_sim_lot]"
    for (size_t i = 0; i < operands; i++) {
        g1_null(g1_operands[i]);
        g1_new(g1_operands[i]);
        g1_rand(g1_operands[i]);
        g2_null(g2_operands[i]);
        g2_new(g2_operands[i]);
        g2_rand(g2_operands[i]);
        bn_new(x[i]);
        bn_null(x[i]);
        bn_rand_mod(x[i], order);
    }

    g1_t temp;
    g1_new(temp);
    g1_null(temp);

    for (size_t i = 0; i < NTESTS; i++) {
        g1_set_infty(temp);
        t[i] = cpucycles();
        for (size_t j = 0; j < operands; j++) {
            g1_mul(g1_operands[j], g1_operands[j], x[j]);
            g1_add(temp, temp, g1_operands[j]);
        }
    }
    printf("[");
    print_results("Results gen param():           ", t, NTESTS);

    for (size_t i = 0; i < NTESTS; i++) {
        g1_set_infty(temp);
        t[i] = cpucycles();
        g1_mul_sim_lot(temp, g1_operands, x, operands);
    }
    print_results("Results gen param():           ", t, NTESTS);

    g2_t temp_g2;
    g2_new(temp_g2);
    g2_null(temp_g2);

    for (size_t i = 0; i < NTESTS; i++) {
        g2_set_infty(temp_g2);
        t[i] = cpucycles();
        for (size_t j = 0; j < operands; j++) {
            g2_mul(g2_operands[j], g2_operands[j], x[j]);
            g2_add(temp_g2, temp_g2, g2_operands[j]);
        }
    }
    print_results("Results gen param():           ", t, NTESTS);

    for (size_t i = 0; i < NTESTS; i++) {
        g2_set_infty(temp_g2);
        t[i] = cpucycles();
        g2_mul_sim_lot(temp_g2, g2_operands, x, operands);
    }
    print_results("Results gen param():           ", t, NTESTS);
    printf("]\n");
    core_clean();
    return 0;
}
