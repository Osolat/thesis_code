//
// Created by benjamin on 2/11/22.
//

#include <cstdio>
#include <iostream>
#include <string>

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

    bn_t x[NTESTS];
    g1_t g1_gen_vars[NTESTS];

    std::cout << "G1 operations" << std::endl;
    for (size_t i = 0; i < NTESTS; i++) {
        bn_null(x[i]);
        bn_new(x[i]);
        bn_rand_mod(x[i], order);
        g1_null(g1_gen_vars[i]);
        g1_new(g1_gen_vars[i]);
        g1_rand(g1_gen_vars[i]);
    }

    g1_t temp_g1;
    g1_null(temp_g1);
    g1_new(temp_g1);

    // g^xi, variable g
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g1_mul(temp_g1, g1_gen_vars[i], x[i]);
    }

    printf("[");
    print_results("Results gen param():           ", t, NTESTS);

    // g^xi, fixed g
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g1_mul_gen(temp_g1, x[i]);
    }
    print_results("Results gen param():           ", t, NTESTS);

    // g_i + g_j
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g1_add(temp_g1, g1_gen_vars[i], g1_gen_vars[(i + 1) % NTESTS]);
    }
    print_results("Results gen param():           ", t, NTESTS);

    for (size_t i = 0; i < NTESTS; i++) {
        bn_free(x[i]);
        g1_free(g1_gen_vars[i]);
    }

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

    // g^x1 * g2^x2
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g1_mul_sim_gen(temp_g1, x_mult_2[i][0], g1_mult_2[i][1], x_mult_2[i][1]);
    }
    print_results("Results gen param():           ", t, NTESTS);

    // g1^x1 * g2^x2
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g1_mul_sim_lot(temp_g1, g1_mult_2[i], x_mult_2[i], 2);
    }
    print_results("Results gen param():           ", t, NTESTS);

    // prod(gi^xi) all varibable. 10 different i
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g1_mul_sim_lot(temp_g1, g1_mult_10[i], x_mult_10[i], 10);
    }
    print_results("Results gen param():           ", t, NTESTS);

    // prod(gi^xi) all varibable. 100 different i

    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g1_mul_sim_lot(temp_g1, g1_mult_100[i], x_mult_100[i], 100);
    }
    print_results("Results gen param():           ", t, NTESTS);

    for (size_t i = 0; i < NTESTS; i++) {
        for (size_t j = 0; j < 2; j++) {
            bn_free(x_mult_2[i][j]);
        }
        for (size_t j = 0; j < 10; j++) {
            bn_free(x_mult_10[i][j]);
        }
        for (size_t j = 0; j < 100; j++) {
            bn_free(x_mult_100[i][j]);
        }
    }

    unsigned char char_arrays_8[NTESTS][8];
    unsigned char char_arrays_16[NTESTS][16];
    unsigned char char_arrays_32[NTESTS][32];
    for (size_t i = 0; i < NTESTS; i++) {
        memcpy(char_arrays_8[i], (void *)memcpy + i, 8);
        memcpy(char_arrays_16[i], (void *)memcpy + i, 16);
        memcpy(char_arrays_32[i], (void *)memcpy + i, 32);
    }

    // Map random 8 byte arrays
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g1_map(temp_g1, char_arrays_8[i], 8);
    }
    print_results("Results gen param():           ", t, NTESTS);

    // Map random 16 byte arrays
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g1_map(temp_g1, char_arrays_16[i], 16);
    }
    print_results("Results gen param():           ", t, NTESTS);

    // Map random 32 byte arrays
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g1_map(temp_g1, char_arrays_32[i], 32);
    }
    print_results("Results gen param():           ", t, NTESTS);

    printf("]\n");

    std::cout << "G2 operations" << std::endl;
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

    // g^xi, variable g
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g2_mul(temp_g2, g2_gen_vars[i], x[i]);
    }

    printf("[");
    print_results("Results gen param():           ", t, NTESTS);

    // g^xi, fixed g
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g2_mul_gen(temp_g2, x[i]);
    }
    print_results("Results gen param():           ", t, NTESTS);

     // g_i + g_j
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g2_add(temp_g2, g2_gen_vars[i], g2_gen_vars[(i + 1) % NTESTS]);
    }
    print_results("Results gen param():           ", t, NTESTS);

    for (size_t i = 0; i < NTESTS; i++) {
        bn_free(x[i]);
        g2_free(g2_gen_vars[i]);
    }

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

    // g^x1 * g2^x2
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g2_mul_sim_gen(temp_g2, x_mult_2[i][0], g2_mult_2[i][1], x_mult_2[i][1]);
    }
    print_results("Results gen param():           ", t, NTESTS);

    // g1^x1 * g2^x2
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g2_mul_sim_lot(temp_g2, g2_mult_2[i], x_mult_2[i], 2);
    }
    print_results("Results gen param():           ", t, NTESTS);

    // prod(gi^xi) all varibable. 10 different i
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g2_mul_sim_lot(temp_g2, g2_mult_10[i], x_mult_10[i], 10);
    }
    print_results("Results gen param():           ", t, NTESTS);

    // prod(gi^xi) all varibable. 100 different i
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g2_mul_sim_lot(temp_g2, g2_mult_100[i], x_mult_100[i], 100);
    }
    print_results("Results gen param():           ", t, NTESTS);

    for (size_t i = 0; i < NTESTS; i++) {
        for (size_t j = 0; j < 2; j++) {
            bn_free(x_mult_2[i][j]);
        }
        for (size_t j = 0; j < 10; j++) {
            bn_free(x_mult_10[i][j]);
        }
        for (size_t j = 0; j < 100; j++) {
            bn_free(x_mult_100[i][j]);
        }
    }

    for (size_t i = 0; i < NTESTS; i++) {
        memcpy(char_arrays_8[i], (void *)memcpy + i, 8);
        memcpy(char_arrays_16[i], (void *)memcpy + i, 16);
        memcpy(char_arrays_32[i], (void *)memcpy + i, 32);
    }

    // Map random 8 byte arrays
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g2_map(temp_g2, char_arrays_8[i], 8);
    }
    print_results("Results gen param():           ", t, NTESTS);

    // Map random 16 byte arrays
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g2_map(temp_g2, char_arrays_16[i], 16);
    }
    print_results("Results gen param():           ", t, NTESTS);

    // Map random 32 byte arrays
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        g2_map(temp_g2, char_arrays_32[i], 32);
    }
    print_results("Results gen param():           ", t, NTESTS);
    printf("]\n");

    std::cout << "Pairing operations" << std::endl;
    // e(g,h) for random g,h
    gt_t gt_temp;
    gt_null(gt_temp);
    gt_new(gt_temp);
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        pc_map(gt_temp, g1_mult_2[i][0], g2_mult_2[i][0]);
    }
    printf("[");
    print_results("Results gen param():           ", t, NTESTS);

    // e(g,h) for random g,h
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        pc_map(gt_temp, g1_mult_2[i][0], g2_mult_2[i][0]);
    }
    print_results("Results gen param():           ", t, NTESTS);

    // Prod(e(g_i, h_i)) for 10 i
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        pc_map_sim(gt_temp, g1_mult_10[i], g2_mult_10[i], 10);
    }
    print_results("Results gen param():           ", t, NTESTS);

    // Prod(e(g_i, h_i)) for 100 i
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        pc_map_sim(gt_temp, g1_mult_100[i], g2_mult_100[i], 100);
    }
    print_results("Results gen param():           ", t, NTESTS);
    printf("]\n");

    std::cout << "GT operations" << std::endl;
    gt_t gt_gen_vars[NTESTS];

    for (size_t i = 0; i < NTESTS; i++) {
        bn_null(x[i]);
        bn_new(x[i]);
        bn_rand_mod(x[i], order);
        gt_null(gt_gen_vars[i]);
        gt_new(gt_gen_vars[i]);
        gt_rand(gt_gen_vars[i]);
    }

    // g^xi, variable g
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        gt_exp(gt_temp, gt_gen_vars[i], x[i]);
    }

    printf("[");
    print_results("Results gen param():           ", t, NTESTS);

    // g^xi, fixed g
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        gt_exp_gen(gt_temp, x[i]);
    }
    print_results("Results gen param():           ", t, NTESTS);

    gt_t gt_mult_2[NTESTS][2];

    for (size_t i = 0; i < NTESTS; i++) {
        for (size_t j = 0; j < 2; j++) {
            bn_null(x_mult_2[i][j]);
            bn_new(x_mult_2[i][j]);
            bn_rand_mod(x_mult_2[i][j], order);
            gt_null(gt_mult_2[i][j]);
            gt_new(gt_mult_2[i][j]);
            gt_rand(gt_mult_2[i][j]);
        }
    }

    // g1^x1 * g2^x2
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        gt_exp_sim(gt_temp, gt_mult_2[i][0], x_mult_2[i][0], gt_mult_2[i][1], x_mult_2[i][1]);
    }

    print_results("Results gen param():           ", t, NTESTS);

    for (size_t i = 0; i < NTESTS; i++) {
        for (size_t j = 0; j < 2; j++) {
            bn_free(x_mult_2[i][j]);
        }
    }
    printf("]\n");
    core_clean();
    return 0;
}
