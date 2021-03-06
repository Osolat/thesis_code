//
// Created by jonas on 2/18/22.
//

//#include "kLin_KP.h"
#include "../lib/k_lin/k_lin_util.h"

#include <iostream>
#include <cstdio>
#include <string>
#include "../bench_defs.h"


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
//unsigned long long resultArray[4];

int main(int argc, char **argv) {
    std::cout << "Benchmarking KP-ABE_UB from K-Lin on attr=" << atoi(argv[1]) << " and k=" << kss <<"\n";
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

    bn_t mul_attributes[test_attr];
    for (int i = 0; i < N_ATTR; ++i) {
        init_null_new_bn_t_var(mul_attributes[i]);
        bn_set_dig(mul_attributes[i], 1);
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
    //float progress = 0.0;
    for (int jo = 0; jo < 1; jo++) {
        //progressBar(100,progress);
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

        bn_t *Av; bn_t *AW; bn_t *AW0; bn_t *AW1;
        bn_t output[kss]; bn_t outputAW[kss * kss]; bn_t outputAW0[kss * kss]; bn_t outputAW1[kss * kss];
        gt_t map_tmp[kss];

        //bn_t one_as_bn;
        //init_null_new_bn_t_var(one_as_bn);
        //bn_set_dig(one_as_bn, 1);

        //TODO fix this, actually dont need it. Just use g1_set_infty?
        //g1_t one_as_g1;
        //init_null_new_g1_t_var(one_as_g1);
        //g1_mul(one_as_g1, group1, one_as_bn);

        //TODO fix these to be optimal
        Av = matrix_mul_vector(output, A1_tmp, msk.v_secret, kss, two_k, two_k, order);
        //AW = matrixG1_mul_matrixBN(outputAW, mpk.A1_mat, msk.W_matrix, kss, two_k, two_k, kss, one_as_g1);
        //AW0 = matrixG1_mul_matrixBN(outputAW0, mpk.A1_mat, msk.W0_matrix, kss, two_k, two_k, kss,one_as_g1);
        //AW1 = matrixG1_mul_matrixBN(outputAW1, mpk.A1_mat, msk.W1_matrix, kss, two_k, two_k, kss, one_as_g1);

        AW = matrixA_mul_matrixW(outputAW, A1_tmp, msk.W_matrix, kss, two_k, two_k, kss, order);
        AW0 = matrixA_mul_matrixW(outputAW0, A1_tmp, msk.W0_matrix, kss, two_k, two_k, kss, order);
        AW1 = matrixA_mul_matrixW(outputAW1, A1_tmp, msk.W1_matrix, kss, two_k, two_k, kss, order);

        for (int k = 0; k < (kss * kss); k++) {
            if (k < kss) {
                //Initialize the gt entries of the e-mapping matrix doing matrix multiplications exponent-wise.
                pp_map_oatep_k12(map_tmp[k], group1, group2);
                gt_exp(mpk.e_mat[k], map_tmp[k], Av[k]);
            }
            g1_mul(mpk.AW_mat[k], group1, AW[k]);
            g1_mul(mpk.AW0_mat[k], group1, AW0[k]);
            g1_mul(mpk.AW1_mat[k], group1, AW1[k]);
        }
        //progress = ((float) (jo+1) / NTESTS);
    }
    //test_stuff(resultArray, 0, t, NTESTS);
    //printf("[");
    //print_results("Results gen param():           ", t, NTESTS);


    /* Key Generation */
    //float progress2 = 0.0;
    struct secret_key_K_Lin_lu sk;
    struct sk_tmp_vectors_lu vj;
    struct node tree_root;

    std::vector <policy_coefficient> res;
    init_secret_key_K_Lin_lu(N_ATTR, &sk);
    init_sk_tmp_vectors_lu(N_ATTR, kss, &vj);
    tree_from_string(and_tree_formula(N_ATTR), &tree_root);

    for (int no = 0; no < NTESTS; no++) {
        //progressBar(100,progress2);
        t[no] = cpucycles();
        bn_t *Wr; bn_t *jW1; bn_t *W0_W1; bn_t *W0_w1_rj; bn_t *v_plus_w;
        bn_t output1[two_k];
        bn_t output3[two_k * kss];
        bn_t output4[two_k * kss];
        bn_t output5[two_k];
        bn_t output1_v_plus_w[two_k];

        for (int i = 0; i < (two_k); ++i) {
            res = std::vector<policy_coefficient>();
            share_secret(&tree_root, msk.v_secret[i], order, res, true);
            for (auto it = res.begin(); it != res.end(); ++it) {
                bn_copy(vj.vj[it->leaf_index - 1].vec_j[i], it->share);
                //g2_mul(sk.sk4[it->leaf_index - 1].sk_four[i], group2, vj.vj[it->leaf_index - 1].vec_j[i]);
            }
        }

        for (auto it = res.begin(); it != res.end(); ++it) {
            for (int k = 0; k < (kss); k++) {
                bn_rand_mod(vj.rj[it->leaf_index - 1].vec_rj[k], order);
                g2_mul(sk.sk13[it->leaf_index - 1].sk_two[k], group2, vj.rj[it->leaf_index - 1].vec_rj[k]);
            }

            jW1 = matrix_mul_scalar(output3, msk.W1_matrix, it->leaf_index - 1, two_k, kss, order);                                                               //h+1 so that we don't multiply with 0
            W0_W1 = matrix_add_matrix(output4, msk.W0_matrix, jW1, two_k, kss, two_k, kss,order);                                     //TODO Maybe not correct with rho(j) !=0 / =0
            W0_w1_rj = matrix_mul_vector(output5, W0_W1, vj.rj[it->leaf_index - 1].vec_rj, two_k, kss, kss, order);

            Wr = matrix_mul_vector(output1, msk.W_matrix, vj.rj[it->leaf_index - 1].vec_rj, two_k, kss, kss, order);
            v_plus_w = vector_add_vector(output1_v_plus_w, vj.vj[it->leaf_index - 1].vec_j, Wr, two_k, two_k, order);

            for (int s = 0; s < (two_k); ++s) {
                g2_mul(sk.sk13[it->leaf_index - 1].sk_three[s], group2, W0_w1_rj[s]);
                g2_mul(sk.sk13[it->leaf_index - 1].sk_one[s], group2, v_plus_w[s]);
            }
        }

        /*
        for (int i = 0; i < (two_k); ++i) {
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
        */
        //progress2 = ((float) (no+1) / NTESTS);
    }
    //test_stuff(resultArray, 1, t, NTESTS);
    printf("[");
    print_results("Results keyGen():           ", t, NTESTS);

    /* Encryption */
    //float progress3 = 0.0;
    //Initialize ciphertext struct
    struct ciphertext_K_Lin_lu CT_A;
    struct tmp_si_lu si;
    init_ciphertext_K_Lin_lu(N_ATTR, kss, &CT_A);
    init_tmp_si_lu(N_ATTR, kss, &si);
    bn_t rnd_s[kss];

    for (int qo = 0; qo < NTESTS; qo++) {
        //progressBar(100,progress3);
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

            for (int p = 0; p < two_k; ++p) {
                if (p < kss) {
                    g1_copy(CT_A.C_23[z].c_2_vec[p], res_Add_val[p]);
                }
                g1_copy(CT_A.C_23[z].c_3_vec[p], ct_3[p]);
            }
        }
        //progress3 = ((float) (qo+1) / NTESTS);
    }
    //test_stuff(resultArray, 2, t, NTESTS);
    print_results("Results encryption():           ", t, NTESTS);

    //TODO start/complete decryption.
    //
    /* Decryption */
    //float progress4 = 0.0;
    //bn_t pack_coef_neg[N_ATTR];
    gt_t map_tmp_1;
    gt_t map_tmp_2;
    gt_t map_tmp_3;
    gt_t map_tmp_4;
    gt_t map_tmp_prod_1;
    gt_t map_tmp_prod_2;
    gt_t map_tmp_prod_3;
    gt_t map_tmp_prod_4;
    gt_t invert_elem_1;
    gt_t de_nom;
    gt_t map_res;

    gt_t exp_val; gt_t exp_val_extra; gt_t prod; gt_t mul_val; gt_t mul_val_extra; gt_t tmp_res;
    init_null_new_gt_t_var(exp_val);
    init_null_new_gt_t_var(exp_val_extra);
    init_null_new_gt_t_var(prod);
    init_null_new_gt_t_var(mul_val);
    init_null_new_gt_t_var(mul_val_extra);
    init_null_new_gt_t_var(tmp_res);
    init_null_new_gt_t_var(map_tmp_1);
    init_null_new_gt_t_var(map_tmp_2);
    init_null_new_gt_t_var(map_tmp_3);
    init_null_new_gt_t_var(map_tmp_4);
    init_null_new_gt_t_var(map_tmp_prod_1);
    init_null_new_gt_t_var(map_tmp_prod_2);
    init_null_new_gt_t_var(map_tmp_prod_3);
    init_null_new_gt_t_var(map_tmp_prod_4);
    init_null_new_gt_t_var(invert_elem_1);
    init_null_new_gt_t_var(de_nom);
    init_null_new_gt_t_var(map_res);

    for (int go = 0; go < NTESTS; go++) {
        //progressBar(100,progress4);
        t[go] = cpucycles();

        gt_set_unity(mul_val);
        gt_set_unity(tmp_res);

        try {
            check_satisfiability(&tree_root, attributes, N_ATTR);
        } catch (struct TreeUnsatisfiableException *e) {
            printf("Fail");
        }

        res = recover_coefficients(&tree_root, attributes, N_ATTR);

        for (auto it5 = res.begin(); it5 != res.end(); ++it5) {
                gt_set_unity(map_tmp_prod_1);
                gt_set_unity(map_tmp_prod_2);
                gt_set_unity(map_tmp_prod_3);
                gt_set_unity(map_tmp_prod_4);

                //TODO refactor this shitty way of doing this. Might not need it when we use all N attributes and with "AND"-tree
                //bn_t neg_coef;
                //init_null_new_bn_t_var(neg_coef);
                //bn_copy(neg_coef, it5->coeff);
                //bn_t_negate(neg_coef, order);
                //bn_copy(pack_coef_neg[it5->leaf_index - 1], neg_coef);

                for (int j = 0; j < ((2 * two_k) + kss); ++j) {
                    if (j < two_k) {
                        pc_map(map_tmp_1, CT_A.C_1[j], sk.sk13[it5->leaf_index - 1].sk_one[j]);
                        gt_mul(map_tmp_prod_1, map_tmp_prod_1, map_tmp_1);
                    } else if (j >= two_k && j < (2 * two_k)) {
                        pc_map(map_tmp_2, CT_A.C_23[it5->leaf_index - 1].c_3_vec[j % two_k], sk.sk13[it5->leaf_index - 1].sk_three[j % two_k]);
                        gt_mul(map_tmp_prod_2, map_tmp_prod_2, map_tmp_2);
                    } else {
                        pc_map(map_tmp_3, CT_A.C_23[it5->leaf_index - 1].c_2_vec[j % (2 * two_k)], sk.sk13[it5->leaf_index - 1].sk_two[j % (2 * two_k)]);
                        gt_mul(map_tmp_prod_3, map_tmp_prod_3, map_tmp_3);
                    }
                }

                //TODO not used right now as we only use policies of AND gates.
                //for (int ik = 0; ik < (two_k); ++ik) {
                //pc_map(map_tmp_4, CT_A.C_1[ik], sk.sk4[it5->leaf_index - 1].sk_four[ik]);
                //gt_mul(map_tmp_prod_4, map_tmp_prod_4, map_tmp_4);
                //}

                gt_mul(de_nom, map_tmp_prod_1, map_tmp_prod_2);
                gt_inv(invert_elem_1, de_nom);
                gt_mul(map_res, invert_elem_1, map_tmp_prod_3);

                //Here we do map_sim = [-sTAv_j]^(wj) where map_sim = [-sTAv_j] comes from the correctness of the K_Lin paper and wj is the coefficients.
                gt_exp(exp_val, map_res, it5->coeff);

                //TODO use for rho(j)=0 in the last mapping of ct1 and ct4
                //gt_exp(exp_val_extra, map_tmp_prod_4, pack_coef_neg[it5->leaf_index - 1]);

                //TODO multiply cases where X_rho(j)=1 and rho(j)=0
                //gt_mul(prod, exp_val, exp_val_extra);

                gt_mul(mul_val, mul_val, exp_val);

        }

        /*
        for (auto it6 = res.begin(); it6 != res.end(); ++it6) {
            if ((it6->leaf_index - 1) == 0) {
                gt_set_unity(map_tmp_prod_4);
                bn_t neg_coef;
                init_null_new_bn_t_var(neg_coef);
                bn_copy(neg_coef, it6->coeff);
                bn_t_negate(neg_coef, order);
                bn_copy(pack_coef_neg[it6->leaf_index - 1], neg_coef);

                //TODO not used right now as we only use policies of AND gates.
                for (int ik = 0; ik < (two_k); ++ik) {
                    pc_map(map_tmp_4, CT_A.C_1[ik], sk.sk4[it6->leaf_index - 1].sk_four[ik]);
                    gt_mul(map_tmp_prod_4, map_tmp_prod_4, map_tmp_4);
                }

                //TODO use for rho(j)=0 in the last mapping of ct1 and ct4
                gt_exp(exp_val_extra, map_tmp_prod_4, pack_coef_neg[it6->leaf_index - 1]);

                //TODO multiply cases where X_rho(j)=1 and rho(j)=0
                //gt_mul(prod, exp_val, exp_val_extra);
                gt_mul(prod, prod, exp_val_extra);
            }
        }
        */



        //Here we complete the product of [-sTAv_j]^(wj)
        gt_mul(tmp_res, tmp_res, mul_val);

        //Printouts for correctness.
        gt_t final_final_res;
        init_null_new_gt_t_var(final_final_res);
        gt_mul(final_final_res, tmp_res, CT_A.C_4_one_val);

        if (!gt_cmp(final_final_res, CT_A.M) == RLC_EQ) {
            printf("Decryption failed!: %d\n", gt_cmp(final_final_res, CT_A.M) == RLC_EQ);
        }

        //progress4 = ((float) (go+1) / NTESTS);
    }
    //test_stuff(resultArray, 3, t, NTESTS);
    //print_result_array(resultArray);

    print_results("Results decryption():           ", t, NTESTS);
    printf("]\n");
    std::cout<<"\n"<<std::endl;
    return 0;
}


