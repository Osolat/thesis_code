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
    std::cout << "[g1 mul, g1 mul fix, g2 mul, g2 mul fix]" << std::endl;
    gt_t m;
    gt_new(m);
    gt_null(m);

    g1_t g1_operands[operands];
    g2_t g2_operands[operands];
    bn_t x[operands];
    for (size_t i = 0; i < operands; i++) {
        g1_null(g1_operands[i]);
        g1_new(g1_operands[i]);
        g1_rand(g1_operands[i]);
        
        g2_null(g2_operands[i]);
        g2_new(g2_operands[i]);
        g2_rand(g2_operands[i]);
        
        bn_null(x[i]);
        bn_new(x[i]);
        bn_rand_mod(x[i], order);
    }

    g1_t pre_g1[operands][RLC_EP_TABLE_MAX];
    g2_t pre_g2[operands][RLC_EP_TABLE_MAX];

    for (size_t i = 0; i < operands; i++) {
        for (size_t j = 0; j < RLC_EP_TABLE_MAX; j++) {
            /* code */
            g1_null(pre_g1[i][j]);
            g1_new(pre_g1[i][j]);
            g2_null(pre_g2[i][j]);
            g2_new(pre_g2[i][j]);
        }
        g1_mul_pre(pre_g1[i], g1_operands[i]);
        g2_mul_pre(pre_g2[i], g2_operands[i]);
    }

    g1_t temp;
    g1_null(temp);
    g1_new(temp);
    g2_t temptwo;
    g2_null(temptwo);
    g2_new(temptwo);

    for (size_t i = 0; i < NTESTS; i++) {
        g1_set_infty(temp);
        t[i] = cpucycles();
        for (size_t j = 0; j < operands; j++) {
            g1_mul(g1_operands[j], g1_operands[j], x[j]);
        }
    }
    printf("[");
    print_results("Results gen param():           ", t, NTESTS);

    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        for (size_t j = 0; j < operands; j++) {
            /* code */
            g1_mul_fix(temp, pre_g1[j], x[j]);
        }
    }
    print_results("Results gen param():           ", t, NTESTS);

    for (size_t i = 0; i < NTESTS; i++) {
        g1_set_infty(temp);
        t[i] = cpucycles();
        for (size_t j = 0; j < operands; j++) {
            g2_mul(g2_operands[j], g2_operands[j], x[j]);
        }
    }
    print_results("Results gen param():           ", t, NTESTS);

    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        for (size_t j = 0; j < operands; j++) {
            /* code */
            g2_mul_fix(temptwo, pre_g2[j], x[j]);
        }
    }
    print_results("Results gen param():           ", t, NTESTS);

    printf("]\n");
    core_clean();
    return 0;
}
