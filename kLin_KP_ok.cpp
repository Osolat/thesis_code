//
// Created by jonas on 3/23/22.
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

unsigned long long t[NTESTS];


int main(int argc, char **argv) {
    std::cout << "Benchmarking KP-ABE from K-Lin\n";


    if (argc == 1) {
        printf("Need to give argument\n");
        return 0;
    }

    int test_attr = atoi(argv[1]);

    for (int iters = 1; iters < (test_attr+1); ++iters) {
        srand(time(NULL));
        std::string keyInput = "";
        std::string encInput = "";
        std::string keyInputWrong = "";

        uint32_t N_ATTR = iters;

        uint32_t *attr_int_list = NULL;
        attr_int_list = (uint32_t *) malloc(sizeof(uint32_t) * iters);

        int d = 1;

        for (int k = 0; k < iters; k++) {
            keyInput = keyInput + "attr" + std::to_string(d);
            encInput = encInput + "attr" + std::to_string(d);

            if (k < iters - 1) {
                keyInput = keyInput + "|";
                encInput = encInput + " and ";
            }
            attr_int_list[k] = d;
            d++;
        }

        keyInputWrong = "attr1|attr2|attr3|attr4";

        bn_t attributes[iters];
        for (int i = 0; i < N_ATTR; ++i) {
            init_null_new_bn_t_var(attributes[i]);
            bn_set_dig(attributes[i], i + 1);
        }

        std::cout << keyInput;
        //std::cout << keyInputWrong;

        struct master_key_k_lin_ok msk;
        struct public_key_k_lin_ok mpk;

        init_master_key_k_lin_ok(N_ATTR, kss, &msk);
        init_public_key_k_lin_ok(N_ATTR, kss, &mpk);

        core_init();

        bn_t order;
        pc_param_set_any();
        pc_param_print();

        g1_get_ord(order);
        //DEBUG UTIL:
        //test_matrix_mul_vector(kss, order);
        //test_vector_trans_mul_matrix(kss, order);
        //test_matrix_mul_matrix(kss, order);
        //test_vector_dot_product(kss, order);

        /* Generate pre-computation tables for g, h */
        g1_t t_pre_g[RLC_EP_TABLE_MAX];
        g2_t t_pre_h[RLC_EP_TABLE_MAX];

        for (int i = 0; i < RLC_EP_TABLE; i++) {
            init_null_new_g1_t_var(t_pre_g[i]);
            init_null_new_g2_t_var(t_pre_h[i]);
        }
        g1_t group1;
        g2_t group2;
        init_null_new_g1_t_var(group1);
        init_null_new_g2_t_var(group2);

        g1_rand(group1);
        g2_rand(group2);
        g1_mul_pre(t_pre_g, group1);
        g2_mul_pre(t_pre_h, group2);

        /* Setup */
        for (int jo = 0; jo < NTESTS; jo++) {
            t[jo] = cpucycles();
            bn_t A_tmp[(kss + 1) * kss];
            //Initializes the v-vector and sets the entries to some random bn_t value modulo the order.
            for (int d = 0; d < ((kss + 1) * kss); ++d) {
                if (d < (kss + 1)) {
                    bn_rand_mod(msk.v_share[d], order);
                }
                //Initializes the bn_t entries of the A-matrix as just random bn_t value modulo the order.
                //Also initializes the g1 entries of the A-matrix by doing matrix multiplications and sets the A_(i,j) to g1^(AW_(i,j)).
                bn_rand_mod(A_tmp[d], order);
                g2_mul_fix(mpk.a_mat[d], t_pre_h, A_tmp[d]);
            }
            //Define dimensions of matrix A and vector v and calculate the matrix-vector product Av.
            int a_rows = kss;
            int a_cols = kss + 1;
            int v_rows = kss + 1;

            bn_t *Av;
            bn_t output[kss];
            Av = matrix_mul_vector(output, A_tmp, msk.v_share, a_rows, a_cols, v_rows, order);
            gt_t map_tmp[kss];

            for (int k = 0; k < kss; k++) {
                //Initialize the gt entries of the e-mapping matrix doing matrix multiplications exponent-wise.
                pp_map_oatep_k12(map_tmp[k], group1, group2);
                gt_exp(mpk.e_mat[k], map_tmp[k], Av[k]);
            }

            //Initializes the "n" W-matrices (master secret key) by setting every entry in these matrices to some random bn value mod the order.
            //In the K-Lin paper the master secret key consists of w_1,...,w_n and describe w_0 = 0.
            for (int j = 0; j < (N_ATTR + 1); j++) {
                for (int m = 0; m < ((kss + 1) * kss); m++) {
                    if (j != 0) {
                        //Sets entries for Wi where i = 1,...., N_att +1 to random bn_t elements
                        bn_rand_mod(msk.atts[j].w[m], order);
                    } else {
                        //Set all entries in W0 to be zero
                        bn_zero(msk.atts[j].w[m]);
                    }
                }
                //Here matrix multiplication is being calculated.

                bn_t one_as_bn;
                init_null_new_bn_t_var(one_as_bn);
                bn_set_dig(one_as_bn, 1);

                g2_t one_as_g2;
                init_null_new_g2_t_var(one_as_g2);
                g2_mul_fix(one_as_g2, t_pre_h, one_as_bn);

                g2_t *AWi;
                g2_t output[kss * kss];
                AWi = matrixG2_mul_matrixBN(output, mpk.a_mat, msk.atts[j].w, kss, (kss + 1), (kss + 1), kss,one_as_g2);

                //Initializes the "n" AW_i (masker public key).
                for (int x = 0; x < (kss * kss); ++x) {
                    g2_copy(mpk.mats[j].w[x], AWi[x]);
                }
            }
        }

        printf("[");
        print_results("Results gen param():           ", t, NTESTS);

        /* Key Generation */
        struct secret_key_K_Lin_ok sk;
        struct sk_tmp_vj_ok vj;
        init_secret_key_K_Lin_ok(N_ATTR, &sk);
        init_sk_tmp_vj_ok(N_ATTR, kss, &vj);

        struct node tree_root;
        tree_from_string(and_tree_formula(N_ATTR), &tree_root);
        std::vector <policy_coefficient> res;

        for (int no = 0; no < NTESTS; no++) {
            t[no] = cpucycles();
            bn_t *Wr;
            bn_t output1[kss + 1];
            int w_rows = (kss + 1);
            int w_cols = kss;
            int r_rows = kss;

            //For all kss+1 secrets in v:
            for (int i = 0; i < (kss + 1); ++i) {
                res = std::vector<policy_coefficient>();
                share_secret(&tree_root, msk.v_share[i], order, res, true);
                for (auto it2 = res.begin(); it2 != res.end(); ++it2) {
                    //Create and set r_j which is a vector of size k of random elements g2 elements, and sets sk_2j = r_j
                    for (int k = 0; k < (kss); k++) {
                        bn_rand_mod(vj.rj[it2->leaf_index - 1].vec_rj[k], order);
                        g1_mul_fix(sk.sk[it2->leaf_index - 1].sk_two[k], t_pre_g, vj.rj[it2->leaf_index - 1].vec_rj[k]);
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
                Wr = matrix_mul_vector(output1, msk.atts[kj + 1].w, vj.rj[kj].vec_rj, w_rows, w_cols, r_rows,order);
                v_plus_w = vector_add_vector(output1_v_plus_w, vj.vj[kj].vec_j, Wr, (kss + 1), (kss + 1), order);

                //Sets sk_1j by adding all vj vectors with the resulting Wr vectors.
                for (int u = 0; u < (kss + 1); ++u) {
                    g1_mul_fix(sk.sk[kj].sk_one[u], t_pre_g, v_plus_w[u]);
                }
            }
        }

        print_results("Results keyGen():           ", t, NTESTS);

        /* Encryption */
        //Initialize ciphertext struct
        struct ciphertext_K_Lin_ok CT_A;
        init_ciphertext_K_Lin_ok(N_ATTR, kss, &CT_A);
        bn_t rnd_s[kss];

        for (int qo = 0; qo < NTESTS; qo++) {
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
            ct_1 = vector_trans_mul_matrix_g2(output, rnd_s, mpk.a_mat, kss, kss + 1, kss);

            //Finishing ct_1 by doing the exponentiation of g.
            for (int t = 0; t < (kss + 1); ++t) {
                g2_copy(CT_A.C_1[t], ct_1[t]);
            }

            //set ct_2i
            //For all N_ATTR + 1 calculate sTAW_i and the reason for doing it over N_ATTR+1 opposed to N_ATTR is because W_0 = 0 and is done to support that the lsss map can be rho(j) = 0.
            for (int a = 0; a < (N_ATTR + 1); ++a) {
                g2_t *ct2_i;
                g2_t output[kss];
                ct2_i = vector_trans_mul_matrix_g2(output, rnd_s, mpk.mats[a].w, kss, kss, kss);

                //Finishing c_2i, by doing the exponentiation of g.
                for (int v = 0; v < kss; ++v) {
                    g2_copy(CT_A.C_2[a].c_2_mat[v], ct2_i[v]);
                }
            }
        }

        print_results("Results encryption():           ", t, NTESTS);
        //TODO start/complete decryption.

        /* Decryption */
        printf("\n");

        //TODO for policies with OR gates this size needs to be changed, according to the K_Lin paper the size would be <=2*N_ATTR
        //List of all the wj coefficients since j = N_ATTR and because we have (kss+1) secrets to be shared from v, the total amount of coefficients becomes N_ATTR * (kss+1).
        bn_t pack_coef[N_ATTR];

        //Two lists used for the pp_map_sim_oatep_k12 operation.
        //Since ct2j and sk2j are vectors of size (kss) and ct1 and sk1j are vectors of size (kss+1) we need to do ((kss+1)+kss) component wise operations.
        g1_t pair_g1[(kss + 1) + kss];
        g2_t pair_g2[(kss + 1) + kss];

        gt_t exp_val;
        gt_t prod;
        init_null_new_gt_t_var(exp_val);
        init_null_new_gt_t_var(prod);

        //Temporary variable supposed to hold intermediate result of the calculations.
        gt_t tmp_res;
        init_null_new_gt_t_var(tmp_res);

        //Initializes the list of coefficients which should yield a size of N_ATTR * (kss+1)
        for (auto it = res.begin(); it != res.end(); ++it) {
            init_null_new_bn_t_var(pack_coef[it->leaf_index - 1]);
        }

        for (int go = 0; go < NTESTS; go++) {
            t[go] = cpucycles();
            gt_t map_sim;
            init_null_new_gt_t_var(map_sim);

            //Sets tmp_mul_list[r] to one so that the multiplication starts out correct.
            fp12_set_dig(prod, 1);

            //Sets tmp_res to one so that the final multiplications starts out correct.
            fp12_set_dig(tmp_res, 1);

            try {
                check_satisfiability(&tree_root, attributes, N_ATTR);
            } catch (struct TreeUnsatisfiableException *e) {
                printf("Fail");
            }

            res = std::vector<policy_coefficient>();
            res = recover_coefficients(&tree_root, attributes, N_ATTR);

            for (auto it3 = res.begin(); it3 != res.end(); ++it3) {
                //Copy all the coefficients to the pack_coef list.
                bn_copy(pack_coef[it3->leaf_index - 1], it3->coeff);
                //Set up the two lists used for the pp_map_sim_oatep_k12 operation
                for (int jk = 0; jk < ((kss + 1) + kss); ++jk) {
                    if (jk < (kss + 1)) {
                        g2_neg(pair_g2[jk], CT_A.C_1[jk]);
                        g1_copy(pair_g1[jk], sk.sk[it3->leaf_index - 1].sk_one[jk]);
                    } else {
                        g2_copy(pair_g2[jk], CT_A.C_2[(it3->leaf_index - 1) + 1].c_2_mat[(jk + 1) % kss]);
                        g1_copy(pair_g1[jk], sk.sk[it3->leaf_index - 1].sk_two[(jk + 1) % kss]);
                    }
                }
                pp_map_sim_oatep_k12(map_sim, pair_g1, pair_g2, ((kss + 1) + kss));
                //Here we do map_sim = [-sTAv_j]^(wj) where map_sim = [-sTAv_j] comes from the correctness of the K_Lin paper and wj is the coefficients.
                gt_exp(exp_val, map_sim, pack_coef[it3->leaf_index - 1]);
                //Here we basically compute the product of [-sTAv_j]^(wj) and saves the result in tmp_mul_list[r]
                gt_mul(prod, prod, exp_val);
            }
            //Here we complete the product of [-sTAv_j]^(wj)
            gt_mul(tmp_res, tmp_res, prod);

            //Printouts for correctness.
            gt_t final_final_res;
            init_null_new_gt_t_var(final_final_res);

            gt_mul(final_final_res, tmp_res, CT_A.C_3_one_val);
            //printf("The final result \n");
            //gt_print(final_final_res);

            //printf("Message M: \n");
            //gt_print(CT_A.M);

            //Uncomment for correctness check;
            assert(gt_cmp(final_final_res, CT_A.M) == RLC_EQ);
            std::cout << "[*] PASSED" << std::endl;
        }

        print_results("Results decryption():           ", t, NTESTS);
        printf("]\n");
    }
    //Test if msk and mpk is initialized correctly:
    //print_msk(&msk, N_ATTR, kss);
    //print_mpk(&mpk, N_ATTR, kss);
    //print_sk(&sk, &vj, N_ATTR, kss);
    //print_ct(&CT_A, N_ATTR, kss);

    return 0;
}

