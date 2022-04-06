//
// Created by jonas on 4/5/22.
//

//
// Created by jonas on 2/18/22.
//

//#include "kLin_KP.h"
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
    std::cout << "Benchmarking KP-ABE_lu from K-Lin\n";
    srand(time(NULL));

    if (argc == 1) {
        printf("Need to give argument\n");
        return 0;
    }

    int test_attr = atoi(argv[1]);
    int two_k = ((2*kss) + 1);

    uint32_t N_ATTR = test_attr;

    bn_t attributes[test_attr];
    for (int i = 0; i < N_ATTR; ++i) {
        init_null_new_bn_t_var(attributes[i]);
        bn_set_dig(attributes[i], i + 1);
    }

    struct master_key_k_lin_lu msk;
    struct public_key_k_lin_lu mpk;

    init_master_key_k_lin_lu(N_ATTR, kss, &msk);
    init_public_key_k_lin_lu(N_ATTR, kss, &mpk);

    core_init();
    bn_t order;
    pc_param_set_any();
    pc_param_print();
    g1_get_ord(order);

    g1_t group1;
    g2_t group2;
    init_null_new_g1_t_var(group1);
    init_null_new_g2_t_var(group2);

    /* Setup */
    for (int jo = 0; jo < 1; jo++) {
        //t[jo] = cpucycles();

        g1_rand(group1);
        g2_rand(group2);

        bn_t A1_tmp[two_k * kss];

        for (int d = 0; d < (two_k * kss); ++d) {
            if (d < (two_k)) {
                bn_rand_mod(msk.v_secret[d], order);
            }
            bn_rand_mod(A1_tmp[d], order);
            bn_rand_mod(msk.W_matrix[d], order);
            bn_rand_mod(msk.W1_matrix[d], order);
            bn_rand_mod(msk.W0_matrix[d], order);
            g1_mul(mpk.A1_mat[d], group1, A1_tmp[d]);
        }

        bn_t *Av; g1_t *AW; g1_t *AW0; g1_t *AW1;
        bn_t output[kss];
        g1_t outputAW[kss * kss];
        g1_t outputAW0[kss * kss];
        g1_t outputAW1[kss * kss];
        gt_t map_tmp[kss];

        bn_t one_as_bn;
        init_null_new_bn_t_var(one_as_bn);
        bn_set_dig(one_as_bn, 1);

        //TODO fix this, actually dont need it. Just use g1_set_infty?
        g1_t one_as_g1;
        init_null_new_g1_t_var(one_as_g1);
        g1_mul(one_as_g1, group1, one_as_bn);

        //TODO fix these to be optimal
        Av = matrix_mul_vector(output, A1_tmp, msk.v_secret, kss, two_k, two_k, order);
        AW = matrixG1_mul_matrixBN(outputAW, mpk.A1_mat, msk.W_matrix, kss, two_k, two_k, kss, one_as_g1);
        AW0 = matrixG1_mul_matrixBN(outputAW0, mpk.A1_mat, msk.W0_matrix, kss, two_k, two_k, kss,one_as_g1);
        AW1 = matrixG1_mul_matrixBN(outputAW1, mpk.A1_mat, msk.W1_matrix, kss, two_k, two_k, kss, one_as_g1);

        for (int k = 0; k < (kss * kss); k++) {
            if (k < kss) {
                //Initialize the gt entries of the e-mapping matrix doing matrix multiplications exponent-wise.
                pp_map_oatep_k12(map_tmp[k], group1, group2);
                gt_exp(mpk.e_mat[k], map_tmp[k], Av[k]);
            }
            g1_copy(mpk.AW_mat[k], AW[k]);
            g1_copy(mpk.AW0_mat[k], AW0[k]);
            g1_copy(mpk.AW1_mat[k], AW1[k]);
        }
    }

    //print_results("Results gen param():           ", t, NTESTS);


    /* Key Generation */
    struct secret_key_K_Lin_lu sk;
    struct sk_tmp_vectors_lu vj;
    struct node tree_root;

    std::vector <policy_coefficient> res;
    init_secret_key_K_Lin_lu(N_ATTR, &sk);
    init_sk_tmp_vectors_lu(N_ATTR, kss, &vj);

    for (int no = 0; no < NTESTS; no++) {
        t[no] = cpucycles();
        free_tree(&tree_root);
        tree_root = node();
        tree_from_string(and_tree_formula(N_ATTR), &tree_root);

        for (int i = 0; i < (two_k); ++i) {
            bn_t *Wr; bn_t *jW1; bn_t *W0_W1; bn_t *W0_w1_rj; bn_t *v_plus_w;
            bn_t output1[two_k];
            bn_t output3[two_k * kss];
            bn_t output4[two_k * kss];
            bn_t output5[two_k];
            bn_t output1_v_plus_w[two_k];

            res = std::vector<policy_coefficient>();
            share_secret(&tree_root, msk.v_secret[i], order, res, true);
            for (auto it = res.begin(); it != res.end(); ++it) {
                //Create and set r_j which is a vector of size k of random elements g2 elements, and sets sk_2j = r_j
                for (int k = 0; k < (kss); k++) {
                    bn_rand_mod(vj.rj[it->leaf_index - 1].vec_rj[k], order);
                    g2_mul(sk.sk13[it->leaf_index - 1].sk_two[k], group2, vj.rj[it->leaf_index - 1].vec_rj[k]);
                }
                //Sets the vj's to contain the j shares for the (kss+1) secrets of v.
                //To clarify each vj is a vector of size (kss+1) and there are a total of j vectors.
                bn_copy(vj.vj[it->leaf_index - 1].vec_j[i], it->share);
                g2_mul(sk.sk4[it->leaf_index - 1].sk_four[i], group2, vj.vj[it->leaf_index - 1].vec_j[i]);                                                                                                         //Correct if rho(j)=0 for all j, however should be empty when all AND gates in policy tree

                jW1 = matrix_mul_scalar(output3, msk.W1_matrix, it->leaf_index - 1, two_k, kss, order);                                                               //h+1 so that we don't multiply with 0
                W0_W1 = matrix_add_matrix(output4, msk.W0_matrix, jW1, two_k, kss, two_k, kss,order);                                     //TODO Maybe not correct with rho(j) !=0 / =0
                W0_w1_rj = matrix_mul_vector(output5, W0_W1, vj.rj[it->leaf_index - 1].vec_rj, two_k, kss, kss, order);

                for (int s = 0; s < (two_k); ++s) {
                    g2_mul(sk.sk13[it->leaf_index - 1].sk_three[s], group2, W0_w1_rj[s]);
                }
            }

            for (auto it2 = res.begin(); it2 != res.end(); ++it2) {
                //Computes W_j * rj by matrix-vector multiplication.
                Wr = matrix_mul_vector(output1, msk.W_matrix, vj.rj[it2->leaf_index - 1].vec_rj, two_k, kss, kss, order);
                v_plus_w = vector_add_vector(output1_v_plus_w, vj.vj[it2->leaf_index - 1].vec_j, Wr, two_k, two_k, order);                                                   //TODO Bug (relic stack smashing) here but only occurs when k=1 and N_Attr=50. Every other values work
                //Sets sk_1j by adding all vj vectors with the resulting Wr vectors.
                for (int u = 0; u < (two_k); ++u) {
                    g2_mul(sk.sk13[it2->leaf_index - 1].sk_one[u], group2, v_plus_w[u]);
                }
            }
        }
    }

    printf("[");
    print_results("Results keyGen():           ", t, NTESTS);

    /* Encryption */
    //Initialize ciphertext struct
    struct ciphertext_K_Lin_lu CT_A;
    struct tmp_si_lu si;
    init_ciphertext_K_Lin_lu(N_ATTR, kss, &CT_A);
    init_tmp_si_lu(N_ATTR, kss, &si);
    bn_t rnd_s[kss];

    for (int qo = 0; qo < NTESTS; qo++) {
        t[qo] = cpucycles();

        gt_t gt_mul_test; gt_t gt_st_test;
        init_null_new_gt_t_var(gt_mul_test);
        init_null_new_gt_t_var(gt_st_test);
        fp12_set_dig(gt_st_test, 1);

        //Set M = 1
        init_null_new_gt_t_var(CT_A.M);
        fp12_set_dig(CT_A.M, 1);

        //Sample random vector s of size k and set c_4
        for (int i = 0; i < kss; ++i) {
            bn_rand_mod(rnd_s[i], order);
            gt_exp(gt_mul_test, mpk.e_mat[i], rnd_s[i]);
            gt_mul(gt_st_test, gt_st_test, gt_mul_test);
        }

        //Multiply gt value of sTAv with the message M to complete ct_4
        gt_mul(CT_A.C_4_one_val, gt_st_test, CT_A.M);


        //set ct_1
        g1_t *ct_1_g1;
        g1_t output_g1[two_k];
        ct_1_g1 = vector_trans_mul_matrix_g1(output_g1, rnd_s, mpk.A1_mat, kss, two_k, kss);

        for (int t = 0; t < (two_k); ++t) {
            g1_copy(CT_A.C_1[t], ct_1_g1[t]);
        }

        for (int z = 0; z < N_ATTR; ++z) {                                                                                                                                              //i=N_ATTR because all attribute is needed to decrypt due the fact it's all AND gates
            for (int x = 0; x < kss; ++x) {
                bn_rand_mod(si.si[z].si_vec[x], order);
            }

            //set ct_3
            g1_t *ct_3;
            g1_t output_ct_3[two_k];

            //Calculate sT*A using vector-matrix multiplication for a transposed vector.
            ct_3 = vector_trans_mul_matrix_g1(output_ct_3, si.si[z].si_vec, mpk.A1_mat, kss, two_k, kss);

            //Finishing ct_1 by doing the exponentiation of g.
            for (int ti = 0; ti < (two_k); ++ti) {
                g1_copy(CT_A.C_23[z].c_3_vec[ti], ct_3[ti]);
            }

            g1_t *sTAW;
            g1_t output[kss];
            sTAW = vector_trans_mul_matrix_g1(output, rnd_s, mpk.AW_mat, kss, kss, kss);

            g1_t *i_A1_W1; g1_t *W0_i_A1_W1; g1_t *res_sT_W0; g1_t *res_Add_val;
            g1_t output_i_A1_W1[kss * kss];
            g1_t output_W0_i_A1_W1[kss * kss];
            g1_t output_res_sT_W0[kss];
            g1_t output_res_Add_val[kss];

            i_A1_W1 = matrix_mul_scalar_g1(output_i_A1_W1, mpk.AW1_mat, (z), kss, kss);
            W0_i_A1_W1 = matrix_add_matrix_g1(output_W0_i_A1_W1, mpk.AW0_mat, i_A1_W1, kss, kss, kss,kss);
            res_sT_W0 = vector_trans_mul_matrix_g1(output_res_sT_W0, si.si[z].si_vec, W0_i_A1_W1, kss, kss, kss);
            res_Add_val = vector_add_vector_g1(output_res_Add_val, sTAW, res_sT_W0, kss, kss);

            for (int p = 0; p < kss; ++p) {
                g1_copy(CT_A.C_23[z].c_2_vec[p], res_Add_val[p]);
            }
        }
    }
    print_results("Results encryption():           ", t, NTESTS);

    //TODO start/complete decryption.
    printf("\n");

    /* Decryption */
    bn_t pack_coef[N_ATTR];
    bn_t pack_coef_neg[N_ATTR];

    g1_t pair_g1[(2 * two_k) + kss];
    g2_t pair_g2[(2 * two_k) + kss];
    g1_t pair_g1_ex[two_k];
    g2_t pair_g2_ex[two_k];

    gt_t exp_val; gt_t exp_val_extra; gt_t prod; gt_t mul_val; gt_t mul_val_extra;
    init_null_new_gt_t_var(exp_val);
    init_null_new_gt_t_var(exp_val_extra);
    init_null_new_gt_t_var(prod);
    init_null_new_gt_t_var(mul_val);
    init_null_new_gt_t_var(mul_val_extra);

    //Temporary variable supposed to hold intermediate result of the calculations.
    gt_t tmp_res;
    init_null_new_gt_t_var(tmp_res);

    //Initializes the list of coefficients which should yield a size of N_ATTR * (kss+1)
    for (auto it4 = res.begin(); it4 != res.end(); ++it4) {
        init_null_new_bn_t_var(pack_coef[it4->leaf_index - 1]);                       //Same as for std.
    }

    for (int go = 0; go < NTESTS; go++) {
        t[go] = cpucycles();

        gt_t map_sim; gt_null(map_sim); gt_new(map_sim);
        gt_t map_sim_2; gt_null(map_sim_2); gt_new(map_sim_2);

        int wj = 0;
        //Sets mul_val to one so that the multiplication starts out correct.
        fp12_set_dig(mul_val, 1);

        //Sets tmp_res to one so that the final multiplications starts out correct.
        fp12_set_dig(tmp_res, 1);

        try {
            check_satisfiability(&tree_root, attributes, N_ATTR);
        } catch (struct TreeUnsatisfiableException *e) {
            printf("Fail");
        }

        res = std::vector<policy_coefficient>();
        res = recover_coefficients(&tree_root, attributes, N_ATTR);

        int j = 0;
        //For all Attributes, set up the two lists used for the pp_map_sim_oatep_k12 operation.
        for (auto it5 = res.begin(); it5 != res.end(); ++it5) {
            //Copy all the coefficients to the pack_coef list.
            bn_copy(pack_coef[wj], it5->coeff);
            //Shitty way of doing this.
            bn_t neg_coef;
            init_null_new_bn_t_var(neg_coef);

            bn_copy(neg_coef, it5->coeff);
            bn_t_negate(neg_coef, order);
            bn_copy(pack_coef_neg[wj], neg_coef);

            //TODO Uncomment these lines below to observe that the coefficients are identical.
            //Set up the two lists used for the pp_map_sim_oatep_k12 operation
            for (int jk = 0; jk < ((2 * two_k) + kss); ++jk) {
                if (jk < two_k) {
                    g1_neg(pair_g1[jk], CT_A.C_1[jk]);
                    g2_copy(pair_g2[jk], sk.sk13[it5->leaf_index - 1].sk_one[jk]);
                } else if (jk >= two_k && jk < (2 * two_k)) {
                    g1_neg(pair_g1[jk], CT_A.C_23[it5->leaf_index - 1].c_3_vec[jk % two_k]);
                    g2_copy(pair_g2[jk], sk.sk13[it5->leaf_index - 1].sk_three[jk % two_k]);
                } else {
                    g1_copy(pair_g1[jk], CT_A.C_23[it5->leaf_index - 1].c_2_vec[jk % (2 * two_k)]);
                    g2_copy(pair_g2[jk], sk.sk13[it5->leaf_index - 1].sk_two[jk % (2 * two_k)]);
                }
            }

            //TODO not used right now as we only use policies of AND gates.
            for (int ik = 0; ik < (two_k); ++ik) {
                g1_copy(pair_g1_ex[ik], CT_A.C_1[ik]);
                g2_copy(pair_g2_ex[ik], sk.sk4[it5->leaf_index - 1].sk_four[ik]);
            }

            //printf("Do this amount of times \n");
            pp_map_sim_oatep_k12(map_sim, pair_g1, pair_g2, ((2 * two_k) + kss));
            pp_map_sim_oatep_k12(map_sim_2, pair_g1_ex, pair_g2_ex, two_k);

            //Here we do map_sim = [-sTAv_j]^(wj) where map_sim = [-sTAv_j] comes from the correctness of the K_Lin paper and wj is the coefficients.
            gt_exp(exp_val, map_sim, pack_coef[wj]);

            //TODO use for rho(j)=0
            gt_exp(exp_val_extra, map_sim_2, pack_coef_neg[wj]);

            gt_mul(prod, exp_val, exp_val_extra);

            gt_mul(mul_val, mul_val, exp_val);

            //TODO multiply cases where X_rho(j)=1 and rho(j)=0
            wj++;
            j++;
        }
        //Here we complete the product of [-sTAv_j]^(wj)
        gt_mul(tmp_res, tmp_res, mul_val);

        //Printouts for correctness.
        gt_t final_final_res;
        init_null_new_gt_t_var(final_final_res);
        gt_mul(final_final_res, tmp_res, CT_A.C_4_one_val);

        //Uncomment for correctness check;
        assert(gt_cmp(final_final_res, CT_A.M) == RLC_EQ);
        std::cout << "[*] PASSED" << std::endl;
    }

    print_results("Results decryption():           ", t, NTESTS);
    printf("]\n");
    return 0;
}


