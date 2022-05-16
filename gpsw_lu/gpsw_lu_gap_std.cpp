//
// Created by benjamin on 2/11/22.
//

#include <cstdio>
#include <string>

#include "../bench_defs.h"

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
    if (argc == 1) {
        printf("Need to give argument\n");
        return 0;
    }

    int test_attr = atoi(argv[1]);
    std::string keyInput = "";
    std::string encInput = "";

    uint32_t N_ATTR = test_attr;

    struct master_key_kp_gpsw_lu_std msk;
    struct public_key_kp_gpsw_lu_std mpk;

    init_master_key_kp_gpsw_lu_std(N_ATTR, &msk);
    init_public_key_kp_gpsw_lu_std(N_ATTR, &mpk);

    core_init();

    bn_t order;
    bn_null(order);
    bn_new(order);
    pc_param_set_any();
    pc_param_print();
    pc_get_ord(order);
    std::cout << "gpsw_lu_gap_std " << N_ATTR << std::endl;


    /* Setup */

    /*Generator of G1*/
    g1_t g;
    g1_null(g);
    g1_new(g);
    g1_get_gen(g);

    g2_t h;
    g2_null(h);
    g2_new(h);
    g2_get_gen(h);

    /*pick y randomly in Z_p*/
    bn_null(msk.y);
    bn_new(msk.y);
    bn_rand_mod(msk.y, order);

    /*Setup PK*/
    g1_mul_gen(mpk.g1, msk.y);
    g2_rand(mpk.g2);

    for (size_t i = 0; i < N_ATTR + 1; i++) {
        g2_rand(mpk.t_values[i]);
    }

    g2_t pre_h1[RLC_EP_TABLE_MAX];
    g1_t pre_g[RLC_EP_TABLE_MAX];
    g1_t pre_g1[RLC_EP_TABLE_MAX];
    for (size_t j = 0; j < RLC_EP_TABLE_MAX; j++) {
        g2_null(pre_h1[j]);
        g2_new(pre_h1[j]);
        g1_null(pre_g[j]);
        g1_new(pre_g[j]);
        g1_null(pre_g1[j]);
        g1_new([pre_g1[j]]);
    }
    g2_mul_pre(pre_h1, mpk.g2);
    g1_mul_pre(pre_g, g);
    g1_mul_pre(pre_g1, mpk.g1);

    /*KeyGeneration*/
    struct secret_key_kp_gpsw_lu_std sk;
    struct node tree_root;
    std::vector<policy_coefficient> res;
    init_secret_key_kp_gpsw_lu_std(N_ATTR, &sk);

    tree_from_string(and_tree_formula(N_ATTR), &tree_root);
    g2_t temp;
    g2_null(temp);
    g2_new(temp);

    bn_t x;
    bn_null(x);
    bn_new(x);

    bn_t r;
    bn_null(r);
    bn_new(r);
    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        /*Secret sharing of y, according to policy tree*/
        res = std::vector<policy_coefficient>();
        share_secret(&tree_root, msk.y, order, res, true);

        for (auto it = res.begin(); it != res.end(); it++) {
            bn_rand_mod(r, order);
            bn_set_dig(x, it->leaf_index);
            t_function_g2_gap(&temp, x, pre_h1, mpk.t_values, N_ATTR, order);
            g2_norm(temp, temp);
            g2_mul_sim(sk.D_values[it->leaf_index - 1], mpk.g2, it->share, temp, r);
            g1_mul_fix(sk.R_values[it->leaf_index - 1], pre_g, r);
        }
    }
    printf("[");
    print_results("Results gen param():           ", t, NTESTS);
    /*Setup precomputation tables for sk*/
    g1_t pre_R_values[N_ATTR][RLC_EP_TABLE_MAX];
    for (size_t i = 0; i < N_ATTR; i++) {
        for (size_t j = 0; j < RLC_EP_TABLE_MAX; j++) {
            g1_null(pre_R_values[i][j]);
            g1_new(pre_R_values[i][j]);
        }
        g1_mul_pre(pre_R_values[i], sk.R_values[i]);
    }

    g2_t pre_D_values[N_ATTR][RLC_EP_TABLE_MAX];
    for (size_t i = 0; i < N_ATTR; i++) {
        for (size_t j = 0; j < RLC_EP_TABLE_MAX; j++) {
            g2_null(pre_D_values[i][j]);
            g2_new(pre_D_values[i][j]);
        }
        g2_mul_pre(pre_D_values[i], sk.D_values[i]);
    }

    /* Encryption */
    // TODO: Fix message construction.
    gt_t message;
    gt_null(message);
    gt_new(message);
    gt_rand(message);
    // gt_print(message);

    bn_t s;
    bn_null(s);
    bn_new(s);

    g1_t g1_prime;
    g1_null(g1_prime);
    g1_new(g1_prime);

    g2_t T;
    g2_null(T);
    g2_new(T);

    struct ciphertext_kp_gpsw_lu_std E;
    init_ciphertext_kp_gpsw_lu_std(test_attr, &E);

    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        bn_rand_mod(s, order);
        g1_mul_fix(g1_prime, pre_g1, s);
        pc_map(E.E_prime, g1_prime, mpk.g2);
        gt_mul(E.E_prime, E.E_prime, message);
        g1_mul_fix(E.E_prime_prime, pre_g, s);

        for (int i = 0; i < N_ATTR; i++) {
            bn_set_dig(x, i + 1);
            t_function_g2_gap(&T, x, pre_h1, mpk.t_values, N_ATTR, order);
            g2_norm(T, T);
            g2_mul(E.E_values[i], T, s);
        }
    }
    print_results("Results gen param():           ", t, NTESTS);

    /*Decryption(E,D) -> message*/
    bn_t attributes[test_attr];
    for (size_t i = 0; i < N_ATTR; i++) {
        bn_null(attributes[i]);
        bn_new(attributes[i]);
        bn_set_dig(attributes[i], i + 1);
    }
    gt_t F_root;
    gt_null(F_root);
    gt_new(F_root);

    gt_t result;
    gt_null(result);
    gt_new(result);

    for (size_t i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();

        try {
            check_satisfiability(&tree_root, attributes, N_ATTR);
        } catch (struct TreeUnsatisfiableException *e) {
            std::cout << e->what() << std::endl;
        }

        res = recover_coefficients(&tree_root, attributes, N_ATTR);


        gt_t upper_map;
        gt_null(upper_map);
        gt_new(upper_map);
        g2_t g2_temp;
        g2_null(g2_temp);
        g2_new(g2_temp);

        g2_t D_vals[res.size()];
        g2_t E_vals[res.size()+1];
        g1_t R_vals[res.size()+1];

        bn_t coeffs[res.size()];

        if (res.size() > g2_pre_sim_switchpoint) {
            for (auto it = res.begin(); it != res.end(); it++) {
                g1_null(R_vals[it->leaf_index - 1]);
                g1_new(R_vals[it->leaf_index - 1]);
                g1_mul_fix(R_vals[it->leaf_index - 1], pre_R_values[it->leaf_index - 1], it->coeff);
                g1_neg(R_vals[it->leaf_index - 1],R_vals[it->leaf_index - 1]);

                g2_null(E_vals[it->leaf_index - 1]);
                g2_new(E_vals[it->leaf_index - 1]);

                g2_null(D_vals[it->leaf_index - 1]);
                g2_new(D_vals[it->leaf_index - 1]);

                bn_null(coeffs[it->leaf_index - 1]);
                bn_new(coeffs[it->leaf_index - 1]);
                bn_copy(coeffs[it->leaf_index - 1], it->coeff);

                g2_copy(D_vals[it->leaf_index - 1], sk.D_values[it->leaf_index - 1]);
                g2_copy(E_vals[it->leaf_index - 1], E.E_values[it->leaf_index - 1]);

            }
            g2_mul_sim_lot(g2_temp, D_vals, coeffs, res.size());
        } else {
            g2_t mul_temp;
            g2_null(mul_temp);
            g2_new(mul_temp);
            g2_set_infty(g2_temp);
            for (auto it = res.begin(); it != res.end(); it++) {
                g1_null(R_vals[it->leaf_index - 1]);
                g1_new(R_vals[it->leaf_index - 1]);
                g1_mul_fix(R_vals[it->leaf_index - 1], pre_R_values[it->leaf_index - 1], it->coeff);
                g1_neg(R_vals[it->leaf_index - 1], R_vals[it->leaf_index - 1]);

                g2_null(E_vals[it->leaf_index - 1]);
                g2_new(E_vals[it->leaf_index - 1]);

                g2_mul_fix(mul_temp, pre_D_values[it->leaf_index - 1], it->coeff);
                g2_add(g2_temp, g2_temp, mul_temp);
                g2_copy(E_vals[it->leaf_index - 1], E.E_values[it->leaf_index - 1]);
            }
        }
        
        g2_copy(E_vals[res.size()], g2_temp);
        g1_copy(R_vals[res.size()], E.E_prime_prime);
        pc_map_sim(F_root, R_vals, E_vals, res.size()+1);
        gt_inv(F_root, F_root);
        gt_mul(result, F_root, E.E_prime);
    }
    print_results("Results gen param():           ", t, NTESTS);
    printf("]\n");

    if (!gt_cmp(message, result) == RLC_EQ) {
        printf("Result of comparison between Message and F_root: %d\n", gt_cmp(message, result) == RLC_EQ);
    }
    return 0;
}
