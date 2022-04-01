//
// Created by jonas on 3/31/22.
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

static void test_stuff(unsigned long long *array, int idx, unsigned long long *t, size_t tlen) {
    std::cout << std::endl;
    size_t i;
    for (i = 0; i < tlen - 1; i++) {
        t[i] = t[i + 1] - t[i];
    }
    array[idx] = average(t, tlen - 1);
}

static void progressBar(int width, float progress) {
    int barWidth = width;

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) + 1 << " %\r";
    std::cout.flush();
}

static void print_result_array(unsigned long long *array) {
    std::cout << std::endl;
    for (int i = 0; i < 4; ++i) {
        printf("%llu \n", array[i]);
    }
}

unsigned long long t[NTESTS];
unsigned long long resultArray[4];

int main(int argc, char **argv) {
    std::cout << "Benchmarking pre_vs_sim mul\n";

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

    g1_t t_pre_g[RLC_EP_TABLE_MAX];
    g2_t t_pre_h[RLC_EP_TABLE_MAX];

    g1_t t_pre_A[test_comp][RLC_EP_TABLE_MAX];
    g2_t t_pre_A_g2[test_comp][RLC_EP_TABLE_MAX];

    for (int i = 0; i < RLC_EP_TABLE_MAX; i++) {
        g1_null(t_pre_g[i]);
        g1_new(t_pre_g[i]);
        g2_null(t_pre_h[i]);
        g2_new(t_pre_h[i]);
    }

    g1_t group1;
    g2_t group2;

    g1_null(group1);
    g1_new(group1);
    g2_null(group2);
    g2_new(group2);

    bn_t rnd_s[test_comp + 1];
    for (int i = 0; i < test_comp; ++i) {
        bn_rand_mod(rnd_s[i], order);
    }

    g1_t mat_g1[test_comp + 1];
    for (int x = 0; x < test_comp; ++x) {
        g1_rand(mat_g1[x]);
    }

    g2_t mat_g2[test_comp + 1];
    for (int x = 0; x < test_comp; ++x) {
        g2_rand(mat_g2[x]);
    }

    for (int i = 0; i < test_comp; ++i) {
        for (int j = 0; j < RLC_EP_TABLE_MAX; ++j) {
            g1_null(t_pre_A[i][j]);
            g1_new(t_pre_A[i][j]);
            g2_null(t_pre_A_g2[i][j]);
            g2_new(t_pre_A_g2[i][j]);
        }
        g1_t rand_elem;
        g2_t rand_elem_g2;

        g1_null(rand_elem);
        g1_new(rand_elem);
        g2_null(rand_elem_g2);
        g2_new(rand_elem_g2);

        g1_rand(rand_elem);
        g2_rand(rand_elem_g2);
        g1_mul_pre(t_pre_A[i], rand_elem);
        g2_mul_pre(t_pre_A_g2[i], rand_elem_g2);
    }

    g1_t tmp_pre;
    g1_t tmp_add_pre;

    g1_null(tmp_pre);
    g1_new(tmp_pre);
    g1_null(tmp_add_pre);
    g1_new(tmp_add_pre);

    g1_set_infty(tmp_add_pre);

    for (int jo = 0; jo < NTESTS; jo++) {
        t[jo] = cpucycles();

        for (int i = 0; i < test_comp; ++i) {
            g1_mul_fix(tmp_pre, t_pre_A[i], rnd_s[i]);
            g1_add(tmp_add_pre, tmp_add_pre, tmp_pre);
        }
    }
    printf("[");
    print_results("Results pre:           ", t, NTESTS);

    g1_t tmp_sim;
    g1_t tmp_add_sim;
    g1_t sim_res;

    g1_null(tmp_sim);
    g1_new(tmp_sim);
    g1_null(tmp_add_sim);
    g1_new(tmp_add_sim);
    g1_null(sim_res);
    g1_new(sim_res);

    g1_set_infty(tmp_add_sim);
    for (int go = 0; go < NTESTS; go++) {
        t[go] = cpucycles();

        for (int i = 0; i < test_comp; ++i) {
            g1_mul_sim(sim_res, mat_g1[i], rnd_s[i + 1], mat_g1[i + 1], rnd_s[i + 1]);
            g1_add(tmp_add_sim, tmp_add_sim, sim_res);
        }
    }
    print_results("Results sim():           ", t, NTESTS);

    g1_null(sim_res);
    g1_new(sim_res);

    for (int jo = 0; jo < NTESTS; jo++) {
        t[jo] = cpucycles();
        g1_mul_sim_lot(sim_res, mat_g1, rnd_s, test_comp + 1);
    }

    print_results("Results sim lot:           ", t, NTESTS);

    g2_t tmp_pre_g2;
    g2_t tmp_add_pre_g2;

    g2_null(tmp_pre_g2);
    g2_new(tmp_pre_g2);
    g2_null(tmp_add_pre_g2);
    g2_new(tmp_add_pre_g2);

    g2_set_infty(tmp_add_pre_g2);

    for (int jo = 0; jo < NTESTS; jo++) {
        t[jo] = cpucycles();

        for (int i = 0; i < test_comp; ++i) {
            g2_mul_fix(tmp_pre_g2, t_pre_A_g2[i], rnd_s[i]);
            g2_add(tmp_add_pre_g2, tmp_add_pre_g2, tmp_pre_g2);
        }
    }

    print_results("Results pre g2:           ", t, NTESTS);

    g2_t tmp_sim_g2;
    g2_t tmp_add_sim_g2;
    g2_t sim_res_g2;

    g2_null(tmp_sim_g2);
    g2_new(tmp_sim_g2);
    g2_null(tmp_add_sim_g2);
    g2_new(tmp_add_sim_g2);
    g2_null(sim_res_g2);
    g2_new(sim_res_g2);

    g2_set_infty(tmp_add_sim_g2);
    for (int go = 0; go < NTESTS; go++) {
        t[go] = cpucycles();

        for (int i = 0; i < test_comp; ++i) {
            g2_mul_sim(sim_res_g2, mat_g2[i], rnd_s[i + 1], mat_g2[i + 1], rnd_s[i + 1]);
            g2_add(tmp_add_sim_g2, tmp_add_sim_g2, sim_res_g2);
        }
    }
    print_results("Results sim g2():           ", t, NTESTS);

    g2_null(sim_res_g2);
    g2_new(sim_res_g2);

    for (int jo = 0; jo < NTESTS; jo++) {
        t[jo] = cpucycles();
        g2_mul_sim_lot(sim_res_g2, mat_g2, rnd_s, test_comp + 1);
    }

    print_results("Results sim lot g2:           ", t, NTESTS);
    printf("]\n");
    std::cout << "\n"
              << std::endl;

    return 0;
}