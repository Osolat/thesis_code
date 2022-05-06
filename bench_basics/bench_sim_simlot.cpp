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
    core_init();

    bn_t order;
    bn_null(order);
    bn_new(order);

    pc_param_set_any();
    pc_param_print();
    pc_get_ord(order);

    g1_t g1_gen_var;
    g2_t g2_gen_var;

    g1_null(g1_gen_var);
    g1_new(g1_gen_var);

    g2_null(g2_gen_var);
    g2_new(g2_gen_var);

    g1_t temp_g1;
    g1_null(temp_g1);
    g1_new(temp_g1); 

    bn_t x[NTESTS];
    g1_t g1_gen_vars[NTESTS];

    bn_t x_mult_2[NTESTS][2];
    g1_t g1_mult_2[NTESTS][2];

    bn_t x_mult_10[NTESTS][10];
    g1_t g1_mult_10[NTESTS][10];

    bn_t x_mult_100[NTESTS][100];
    g1_t g1_mult_100[NTESTS][100];

    for (size_t i = 0; i < NTESTS; i++) {
        for (size_t j = 0; j < 2; j++) {
            bn_null(x_mult_2[i][j]);
            bn_new(x_mult_2[i][j]);
            bn_rand_mod(x_mult_2[i][j], order);
            g1_null(g1_mult_2[i][j]);
            g1_new(g1_mult_2[i][j]);
            g1_rand(g1_mult_2[i][j]);
        }
        for (size_t j = 0; j < 10; j++) {
            bn_null(x_mult_10[i][j]);
            bn_new(x_mult_10[i][j]);
            bn_rand_mod(x_mult_10[i][j], order);
            g1_null(g1_mult_10[i][j]);
            g1_new(g1_mult_10[i][j]);
            g1_rand(g1_mult_10[i][j]);
        }
        for (size_t j = 0; j < 100; j++) {
            bn_null(x_mult_100[i][j]);
            bn_new(x_mult_100[i][j]);
            bn_rand_mod(x_mult_100[i][j], order);
            g1_null(g1_mult_100[i][j]);
            g1_new(g1_mult_100[i][j]);
            g1_rand(g1_mult_100[i][j]);
        }
    }

    // g1^x1 * g2^x2
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g1_mul_sim(temp_g1, g1_mult_2[i][0], x_mult_2[i][0], g1_mult_2[i][1], x_mult_2[i][1]);
    }
    print_results("Results gen param():           ", t, NTESTS);

    // g1^x1 * g2^x2
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g1_mul_sim_lot(temp_g1, g1_mult_2[i], x_mult_2[i], 2);
    }
    print_results("Results gen param():           ", t, NTESTS);

    printf("]\n");

    std::cout << "G2 operations" << std::endl;
    std::cout << "[mul, mul_gen, add, mul_sim, mul_sim_gen, mul_sim_lot(2), mul_sim_lot(10), mul_sim_lot(100), hash(32)]" << std::endl;

    // G2 tests
    g2_t g2_gen_vars[NTESTS];

    for (size_t i = 0; i < NTESTS; i++) {
        bn_null(x[i]);
        bn_new(x[i]);
        bn_rand_mod(x[i], order);
        g2_null(g2_gen_vars[i]);
        g2_new(g2_gen_vars[i]);
        g2_rand(g2_gen_vars[i]);
    }

    g2_t temp_g2;
    g2_null(temp_g2);
    g2_new(temp_g2);


    g2_t g2_mult_2[NTESTS][2];

    g2_t g2_mult_10[NTESTS][10];

    g2_t g2_mult_100[NTESTS][100];

    for (size_t i = 0; i < NTESTS; i++) {
        for (size_t j = 0; j < 2; j++) {
            bn_null(x_mult_2[i][j]);
            bn_new(x_mult_2[i][j]);
            bn_rand_mod(x_mult_2[i][j], order);
            g2_null(g2_mult_2[i][j]);
            g2_new(g2_mult_2[i][j]);
            g2_rand(g2_mult_2[i][j]);
        }
        for (size_t j = 0; j < 10; j++) {
            bn_null(x_mult_10[i][j]);
            bn_new(x_mult_10[i][j]);
            bn_rand_mod(x_mult_10[i][j], order);
            g2_null(g2_mult_10[i][j]);
            g2_new(g2_mult_10[i][j]);
            g2_rand(g2_mult_10[i][j]);
        }
        for (size_t j = 0; j < 100; j++) {
            bn_null(x_mult_100[i][j]);
            bn_new(x_mult_100[i][j]);
            bn_rand_mod(x_mult_100[i][j], order);
            g2_null(g2_mult_100[i][j]);
            g2_new(g2_mult_100[i][j]);
            g2_rand(g2_mult_100[i][j]);
        }
    }

    // g1^x1 * g2^x2
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g2_mul_sim(temp_g2, g2_mult_2[i][0], x_mult_2[i][0], g2_mult_2[i][1], x_mult_2[i][1]);
    }
    print_results("Results gen param():           ", t, NTESTS);

    // g1^x1 * g2^x2
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g2_mul_sim_lot(temp_g2, g2_mult_2[i], x_mult_2[i], 2);
    }
    print_results("Results gen param():           ", t, NTESTS);
    printf("]\n");
    core_clean();
    return 0;
}
