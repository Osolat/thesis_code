//
// Created by jonas on 3/31/22.
//

#include "lib/k_lin/k_lin_util.h"

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
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
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
    std::cout << "Benchmarking KP-ABE from K-Lin_GAP_OD on attr=" << atoi(argv[1]) << " and k=" << kss <<"\n";

    if (argc == 1) {
        printf("Need to give argument\n");
        return 0;
    }

    int test_attr = atoi(argv[1]);

    srand(time(NULL));
    uint32_t N_ATTR = test_attr;

    bn_t attributes[test_attr];
    for (int i = 0; i < N_ATTR; ++i) {
        init_null_new_bn_t_var(attributes[i]);
        bn_set_dig(attributes[i], i + 1);
    }

    struct master_key_k_lin_od msk;
    struct public_key_k_lin_od mpk;

    init_master_key_k_lin_od(N_ATTR, kss, &msk);
    init_public_key_k_lin_od(N_ATTR, kss, &mpk);

    core_init();
    bn_t order;
    pc_param_set_any();
    pc_param_print();
    g1_get_ord(order);

    /* Generate pre-computation tables for g, h */
    g2_t t_pre_A[(kss + 1) * kss][RLC_EP_TABLE_MAX];
    g1_t t_pre_AW[N_ATTR + 1][kss * kss][RLC_EP_TABLE_MAX];

    g1_t group1;
    g2_t group2;
    init_null_new_g1_t_var(group1);
    init_null_new_g2_t_var(group2);
    //printf("here\n");
    /* Setup */
    //float progress = 0.0;
    for (int jo = 0; jo < 1; jo++) {
        //progressBar(100, progress);
        //t[jo] = cpucycles();

        g1_get_gen(group1);
        g2_get_gen(group2);

        bn_t A_tmp[(kss + 1) * kss];
        //Initializes the v-vector and sets the entries to some random bn_t value modulo the order.
        for (int d = 0; d < ((kss + 1) * kss); ++d) {
            if (d < (kss + 1)) {
                bn_rand_mod(msk.v_share[d], order);
            }
            //Initializes the bn_t entries of the A-matrix as just random bn_t value modulo the order.
            //Also initializes the g1 entries of the A-matrix by doing matrix multiplications and sets the A_(i,j) to g1^(AW_(i,j)).
            bn_rand_mod(A_tmp[d], order);
            g2_mul_gen(mpk.a_mat[d], A_tmp[d]);

            for (int j = 0; j < RLC_EP_TABLE_MAX; ++j) {
                init_null_new_g2_t_var(t_pre_A[d][j]);
            }
            g2_mul_pre(t_pre_A[d], mpk.a_mat[d]);
        }

        bn_t *Av;
        bn_t output[kss];
        Av = matrix_mul_vector(output, A_tmp, msk.v_share, kss, kss + 1, kss + 1, order);
        gt_t map_tmp[kss];

        //Initializes the "n" W-matrices (master secret key) by setting every entry in these matrices to some random bn value mod the order.
        //In the K-Lin paper the master secret key consists of w_1,...,w_n and describe w_0 = 0.
        for (int j = 0; j < (N_ATTR + 1); j++) {
            for (int m = 0; m < ((kss + 1) * kss); m++) {
                if (j != 0) {
                    //Sets entries for Wi where i = 1,...., N_att +1 to random bn_t elements
                    bn_rand_mod(msk.atts[j].w[m], order);
                } else if (m < kss) {
                    pp_map_oatep_k12(map_tmp[m], group1, group2);
                    gt_exp(mpk.e_mat[m], map_tmp[m], Av[m]);
                } else {
                    //Set all entries in W0 to be zero
                    bn_zero(msk.atts[j].w[m]);
                }
            }
            //Here matrix multiplication is being calculated.
            bn_t *AWi;
            bn_t output[kss * kss];
            AWi = matrixA_mul_matrixW(output, A_tmp, msk.atts[j].w, kss, (kss + 1), (kss + 1), kss, order);

            //Initializes the "n" AW_i (masker public key).
            for (int x = 0; x < (kss * kss); ++x) {
                g1_mul_gen(mpk.mats[j].w[x], AWi[x]);
                for (int d = 0; d < RLC_EP_TABLE_MAX; ++d) {
                    init_null_new_g1_t_var(t_pre_AW[j][x][d]);
                }
                g1_mul_pre(t_pre_AW[j][x], mpk.mats[j].w[x]);
            }
        }
        //progress = ((float) (jo+1) / NTESTS);
    }
    //std::cout<<"After\n"<<std::endl;
    //test_stuff(resultArray, 0, t, NTESTS);

    //printf("[");
    //print_results("Results gen param():           ", t, NTESTS);

    /* Key Generation */
    //float progress2 = 0.0;

    struct secret_key_K_Lin_od sk;
    struct sk_tmp_vj_od vj;

    struct node tree_root;
    std::vector <policy_coefficient> res;

    init_secret_key_K_Lin_od(N_ATTR, &sk);
    init_sk_tmp_vj_od(N_ATTR, kss, &vj);
    tree_from_string(and_tree_formula(N_ATTR), &tree_root);

    for (int no = 0; no < NTESTS; no++) {
        //progressBar(100,progress2);

        t[no] = cpucycles();
        bn_t *Wr;
        bn_t output1[kss + 1];

        //For all kss+1 secrets in v:
        for (int i = 0; i < (kss + 1); ++i) {
            res = std::vector<policy_coefficient>();
            share_secret(&tree_root, msk.v_share[i], order, res, true);
            for (auto it2 = res.begin(); it2 != res.end(); ++it2) {
                //Create and set r_j which is a vector of size k of random elements g2 elements, and sets sk_2j = r_j
                for (int k = 0; k < (kss); k++) {
                    bn_rand_mod(vj.rj[it2->leaf_index - 1].vec_rj[k], order);
                    g2_mul_gen(sk.sk[it2->leaf_index - 1].sk_two[k], vj.rj[it2->leaf_index - 1].vec_rj[k]);
                }
                //Sets the vj's to contain the j shares for the (kss+1) secrets of v.
                //To clarify each vj is a vector of size (kss+1) and there are a total of j vectors.
                bn_copy(vj.vj[it2->leaf_index - 1].vec_j[i], it2->share);
            }
        }

        bn_t *v_plus_w;
        bn_t output1_v_plus_w[kss + 1];

        for (int kj = 0; kj < N_ATTR; kj++) {
            //Computes W_j * rj by matrix-vector multiplication.
            Wr = matrix_mul_vector(output1, msk.atts[kj + 1].w, vj.rj[kj].vec_rj, (kss + 1), kss, kss, order);
            v_plus_w = vector_add_vector(output1_v_plus_w, vj.vj[kj].vec_j, Wr, (kss + 1), (kss + 1), order);

            //Sets sk_1j by adding all vj vectors with the resulting Wr vectors.
            for (int u = 0; u < (kss + 1); ++u) {
                g1_mul_gen(sk.sk[kj].sk_one[u], v_plus_w[u]);
            }
        }
        //progress2 = ((float) (no+1) / NTESTS);
    }
    //test_stuff(resultArray, 1, t, NTESTS);
    printf("[");
    print_results("Results keyGen():           ", t, NTESTS);

    /* Encryption */
    //Initialize ciphertext struct
    //float progress3 = 0.0;
    struct ciphertext_K_Lin_od CT_A;
    init_ciphertext_K_Lin_od(N_ATTR, kss, &CT_A);
    bn_t rnd_s[kss];

    for (int qo = 0; qo < NTESTS; qo++) {
        //progressBar(100, progress3);

        t[qo] = cpucycles();
        gt_t gt_mul_test;
        gt_t gt_st_test;
        init_null_new_gt_t_var(gt_mul_test);
        init_null_new_gt_t_var(gt_st_test);
        fp12_set_dig(gt_st_test, 1);

        //Set M = 1
        init_null_new_gt_t_var(CT_A.M);
        fp12_set_dig(CT_A.M, 1);

        //Sample random vector s of size k and set c_3
        for (int i = 0; i < kss; ++i) {
            bn_rand_mod(rnd_s[i], order);
            gt_exp(gt_mul_test, mpk.e_mat[i], rnd_s[i]);
            gt_mul(gt_st_test, gt_st_test, gt_mul_test);
        }

        //Multiply gt value of sTAv with the message M to complete ct_3
        gt_mul(CT_A.C_3_one_val, gt_st_test, CT_A.M);

        //set ct_1
        g2_t *ct_1;
        g2_t output[kss + 1];

        //Calculate sT*A using vector-matrix multiplication for a transposed vector.
        ct_1 = vector_trans_mul_matrix_g2_pre(output, rnd_s, t_pre_A, kss, kss + 1, kss);
        //Finishing ct_1 by doing the exponentiation of g.

        //set ct_2i
        //For all N_ATTR + 1 calculate sTAW_i and the reason for doing it over N_ATTR+1 opposed to N_ATTR is because W_0 = 0 and is done to support that the lsss map can be rho(j) = 0.
        for (int a = 0; a < (N_ATTR + 1); ++a) {
            g1_t *ct2_i;
            g1_t output[kss];
            ct2_i = vector_trans_mul_matrix_g1_pre(output, rnd_s, t_pre_AW[a], kss, kss, kss);

            //Finishing c_2i, by doing the exponentiation of g.
            for (int v = 0; v < (kss + 1); ++v) {
                if (v < kss){
                    g1_copy(CT_A.C_2[a].c_2_mat[v], ct2_i[v]);
                }
                g2_copy(CT_A.C_1[v], ct_1[v]);
            }
        }
        //progress3 = ((float) (qo+1) / NTESTS);
    }
    //test_stuff(resultArray, 2, t, NTESTS);

    print_results("Results encryption():           ", t, NTESTS);

    /* Decryption */
    //float progress4 = 0.0;
    bn_t pack_coef[N_ATTR];
    g1_t g1_list_lol[(N_ATTR * kss) + (kss + 1)];
    g2_t g2_list_lol[(N_ATTR * kss) + (kss + 1)];
    g1_t K1_prod[kss + 1];
    g1_t sk1_tmp[N_ATTR];
    gt_t test_res;

    for (int hg = 0; hg < kss + 1; ++hg) {
        init_null_new_g1_t_var(K1_prod[hg]);
        init_null_new_g1_t_var(sk1_tmp[hg]);
    }
    init_null_new_gt_t_var(test_res);

    gt_t map_sim_test_1;
    init_null_new_gt_t_var(map_sim_test_1);

    for (int go = 0; go < NTESTS; go++) {
        //progressBar(100, progress4);
        t[go] = cpucycles();
        for (auto it4 = res.begin(); it4 != res.end(); ++it4) {
            init_null_new_bn_t_var(pack_coef[it4->leaf_index - 1]);                       //Same as for std.
        }

        try {
            check_satisfiability(&tree_root, attributes, N_ATTR);
        } catch (struct TreeUnsatisfiableException *e) {
            printf("Fail");
        }
        res = std::vector<policy_coefficient>();
        res = recover_coefficients(&tree_root, attributes, N_ATTR);

        int lel = 0;
        int kul = 0;
        for (int po = 0; po < kss + 1; ++po) {
            for (auto it5 = res.begin(); it5 != res.end(); ++it5) {
                if (po == 0) {
                    bn_copy(pack_coef[it5->leaf_index - 1], it5->coeff);
                    for (int jk2 = 0; jk2 < kss; ++jk2) {
                        g1_mul(g1_list_lol[lel], CT_A.C_2[it5->leaf_index].c_2_mat[jk2], pack_coef[it5->leaf_index - 1]);
                        g2_copy(g2_list_lol[lel], sk.sk[it5->leaf_index - 1].sk_two[jk2]);
                        lel++;
                    }
                }
                g1_copy(sk1_tmp[it5->leaf_index - 1], sk.sk[it5->leaf_index - 1].sk_one[po]);
            }
            kul = lel + po;
            g1_mul_sim_lot(K1_prod[po], sk1_tmp, pack_coef, N_ATTR);
            g1_neg(g1_list_lol[kul], K1_prod[po]);
            g2_copy(g2_list_lol[kul], CT_A.C_1[po]);
        }

        pp_map_sim_oatep_k12(map_sim_test_1, g1_list_lol, g2_list_lol, (N_ATTR*kss)+(kss + 1));
        gt_mul(test_res, map_sim_test_1, CT_A.C_3_one_val);

        //Uncomment for correctness check;
        //assert(gt_cmp(test_res, CT_A.M) == RLC_EQ);
        //std::cout << "[*] PASSED" << std::endl;
        //progress4 = ((float) (go+1) / NTESTS);
    }

    //test_stuff(resultArray, 3, t, NTESTS);

    print_results("Results decryption():           ", t, NTESTS);
    printf("]\n");

    std::cout << "\n" << std::endl;
    return 0;
}


