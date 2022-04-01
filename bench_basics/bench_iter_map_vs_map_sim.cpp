//
// Created by jonas on 4/1/22.
//

#include <iostream>
#include <cstdio>
#include <string>
#include "bench_defs.h"


long long cpucycles(void) {
    unsigned long long result;
    asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
    : "=a" (result)::"%rdx");
    return result;
}

static int cmp_llu(const void *a, const void *b) {
    if (*(unsigned long long *) a < *(unsigned long long *) b) return -1;
    if (*(unsigned long long *) a > *(unsigned long long *) b) return 1;
    return 0;
}

static unsigned long long median(unsigned long long *l, size_t llen) {
    qsort(l, llen, sizeof(unsigned long long), cmp_llu);

    if (llen % 2) return l[llen / 2];
    else return (l[llen / 2 - 1] + l[llen / 2]) / 2;
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
    std::cout << "Benchmarking iter_map_vs_map_sim\n";


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

    gt_t iter_res;
    gt_t iter_mul_tmp;
    gt_t copy_sim_res;

    gt_null(iter_res);
    gt_new(iter_res);
    gt_null(iter_mul_tmp);
    gt_new(iter_mul_tmp);
    gt_null(copy_sim_res);
    gt_new(copy_sim_res);

    g1_t g1_ops[test_comp];
    g2_t g2_ops[test_comp];

    g1_t g1_ops_copy[test_comp];
    g2_t g2_ops_copy[test_comp];

    for (int i = 0; i < test_comp; ++i) {

        g1_null(g1_ops[i]);
        g1_new(g1_ops[i]);
        g2_null(g2_ops[i]);
        g2_new(g2_ops[i]);

        g1_null(g1_ops_copy[i]);
        g1_new(g1_ops_copy[i]);
        g2_null(g2_ops_copy[i]);
        g2_new(g2_ops_copy[i]);

        g1_rand(g1_ops[i]);
        g2_rand(g2_ops[i]);
    }

    for (int jo = 0; jo < NTESTS; jo++) {
        t[jo] = cpucycles();
        gt_set_unity(iter_mul_tmp);
        for (int x = 0; x < test_comp; ++x) {
            pc_map(iter_res, g1_ops[x], g2_ops[x]);
            gt_mul(iter_mul_tmp, iter_mul_tmp, iter_res);
        }


    }
    printf("[");
    print_results("Results iter_map:           ", t, NTESTS);



    for (int go = 0; go < NTESTS; go++) {
        t[go] = cpucycles();
        for (int x = 0; x < test_comp; ++x) {
            g1_copy(g1_ops_copy[x], g1_ops[x]);
            g2_copy(g2_ops_copy[x], g2_ops[x]);
        }
        pc_map_sim(copy_sim_res, g1_ops_copy, g2_ops_copy, test_comp);
    }
    print_results("Results sim():           ", t, NTESTS);
    printf("]\n");
    std::cout<<"\n"<<std::endl;

    return 0;
}

