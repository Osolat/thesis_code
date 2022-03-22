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

    std::string keyInput = "";
    std::string encInput = "";
    std::string keyInputWrong = "";

    uint32_t N_ATTR = test_attr;

    uint32_t *attr_int_list = NULL;
    attr_int_list = (uint32_t *) malloc(sizeof(uint32_t) * test_attr);

    int d = 1;

    for (int k = 0; k < test_attr; k++) {
        keyInput = keyInput + "attr" + std::to_string(d);
        encInput = encInput + "attr" + std::to_string(d);

        if (k < test_attr - 1) {
            keyInput = keyInput + "|";
            encInput = encInput + " and ";
        }
        attr_int_list[k] = d;
        d++;
    }

    keyInputWrong = "attr1|attr2|attr3|attr4";

    std::cout << keyInput;

    struct master_key_k_lin_lu msk;
    struct public_key_k_lin_lu mpk;

    init_master_key_k_lin_lu(N_ATTR, kss, &msk);
    init_public_key_k_lin_lu(N_ATTR, kss, &mpk);

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
    //test_matrix_mul_scalar(kss, order);
    //test_matrix_add_matrix(kss, order);
    //test_vector_add_vector(kss, order);

    std::unique_ptr <L_OpenABEFunctionInput> keyFuncInput = nullptr;
    keyFuncInput = L_createAttributeList(keyInput);

    if (keyFuncInput == nullptr) {
        printf("Invalid attribute key input\n");
        return -1;
    }

    L_OpenABEAttributeList *attrList = nullptr;

    if ((attrList = dynamic_cast<L_OpenABEAttributeList *>(keyFuncInput.get())) == nullptr) {
        printf("Error in attribute list\n");
        exit(-1);
    }


    std::unique_ptr <L_OpenABEFunctionInput> keyFuncInputLOL = nullptr;
    keyFuncInputLOL = L_createAttributeList(keyInputWrong);

    if (keyFuncInputLOL == nullptr) {
        printf("Invalid attribute key input\n");
        return -1;
    }

    L_OpenABEAttributeList *attrListLOL = nullptr;

    if ((attrListLOL = dynamic_cast<L_OpenABEAttributeList *>(keyFuncInputLOL.get())) == nullptr) {
        printf("Error in attribute list\n");
        exit(-1);
    }


    /* Generate pre-computation tables for g, h */
    g1_t t_pre_g[RLC_EP_TABLE_MAX];
    g2_t t_pre_h[RLC_EP_TABLE_MAX];

    for (int i = 0; i < RLC_EP_TABLE; i++) {
        g1_new(t_pre_g[i]);
        g2_new(t_pre_h[i]);
    }
    g1_t group1; g1_null(group1); g1_new(group1);
    g2_t group2; g2_null(group2); g2_new(group2);

    g1_rand(group1);
    g2_rand(group2);
    g1_mul_pre(t_pre_g, group1);
    g2_mul_pre(t_pre_h, group2);

    int two_k = ((2*kss) + 1);

    /* Setup */
    for (int jo = 0; jo < NTESTS; jo++) {
        t[jo] = cpucycles();

        bn_t A1_tmp[two_k * kss];


        for (int d = 0; d < (two_k * kss); ++d) {
            if (d < (two_k)) {
                bn_rand_mod(msk.v_secret[d], order);
            }
            bn_rand_mod(A1_tmp[d], order);
            bn_rand_mod(msk.W_matrix[d], order);
            bn_rand_mod(msk.W1_matrix[d], order);
            bn_rand_mod(msk.W0_matrix[d], order);
            g1_mul_fix(mpk.A1_mat[d], t_pre_g, A1_tmp[d]);
        }

        bn_t *Av; g1_t *AW; g1_t *AW0; g1_t *AW1;
        bn_t output[kss];
        g1_t outputAW[kss * kss];
        g1_t outputAW0[kss * kss];
        g1_t outputAW1[kss * kss];
        gt_t map_tmp[kss];

        bn_t one_as_bn;
        bn_null(one_as_bn);
        bn_new(one_as_bn);
        bn_set_dig(one_as_bn, 1);
        g1_t one_as_g1;
        g1_null(one_as_g1);
        g1_new(one_as_g1);
        g1_mul_fix(one_as_g1, t_pre_g, one_as_bn);

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

    printf("[");
    print_results("Results gen param():           ", t, NTESTS);


    /* Key Generation */
    std::unique_ptr <L_OpenABEFunctionInput> funcInput = nullptr;

    struct secret_key_K_Lin_lu sk;
    struct sk_tmp_vectors_lu vj;
    init_secret_key_K_Lin_lu(N_ATTR, &sk);
    init_sk_tmp_vectors_lu(N_ATTR, kss, &vj);


    //TODO this way of setting up OpenAbe lsss seems to work as expected but maybe I am wrong. Check correctness again.
    //TODO for policies with OR gates this size needs to be changed, according to the K_Lin paper the size would be <=2*N_ATTR
    L_OpenABELSSS *lsss_2 = new L_OpenABELSSS(1);
    lsss_2 = (L_OpenABELSSS *) malloc(((N_ATTR+1) * (two_k)) * sizeof(L_OpenABELSSS));
    L_OpenABELSSSRowMap *lsss_row = new L_OpenABELSSSRowMap;
    lsss_row = (L_OpenABELSSSRowMap *) malloc(((N_ATTR+1) * (two_k)) * sizeof(L_OpenABELSSSRowMap) * sizeof(L_OpenABELSSS));

    L_OpenABEPolicy *policy;

    funcInput = L_createPolicyTree(encInput);
    if (!funcInput) {
        printf("Create policy error in encryption\n");
        return -1;
    }

    policy = dynamic_cast< L_OpenABEPolicy *>(funcInput.get());

    if (policy == nullptr) {
        printf("Error in input policy\n");
        return -1;
    }

    for (int no = 0; no < NTESTS; no++) {
        t[no] = cpucycles();


        for (int i = 0; i < (two_k); ++i) {
            L_ZP s_aux;
            s_aux.isOrderSet = true;
            bn_copy(s_aux.order, order);
            bn_copy(s_aux.m_ZP, msk.v_secret[i]);

            //Obtain the j shares for each secret in v. Here we j = N_ATT because the policy consists of only AND gates and so all nodes of the policy tree adds exactly a single share.
            lsss_2[i].l_shareSecret(policy, s_aux);
            lsss_row[i] = lsss_2[i].l_getRows();

            bn_t *Wr; bn_t *jW1; bn_t *W0_W1; bn_t *W0_w1_rj; bn_t *v_plus_w;

            bn_t output1[two_k];
            bn_t output3[two_k * kss];
            bn_t output4[two_k * kss];
            bn_t output5[two_k];
            bn_t output1_v_plus_w[two_k];

            int h = 0;
            for (auto it = lsss_row[i].begin(); it != lsss_row[i].end(); ++it) {
                //Create and set r_j which is a vector of size k of random elements g2 elements, and sets sk_2j = r_j
                for (int k = 0; k < (kss); k++) {
                    bn_rand_mod(vj.rj[h].vec_rj[k], order);
                    g2_mul_fix(sk.sk13[h].sk_two[k], t_pre_h, vj.rj[h].vec_rj[k]);
                }
                //Sets the vj's to contain the j shares for the (kss+1) secrets of v.
                //To clarify each vj is a vector of size (kss+1) and there are a total of j vectors.
                bn_copy(vj.vj[h].vec_j[i], it->second.element().m_ZP);
                g2_mul_fix(sk.sk4[h].sk_four[i], t_pre_h, vj.vj[h].vec_j[i]);                                                                                                         //Correct if rho(j)=0 for all j, however should be empty when all AND gates in policy tree

                jW1 = matrix_mul_scalar(output3, msk.W1_matrix, (h), two_k, kss, order);                                                               //h+1 so that we don't multiply with 0
                W0_W1 = matrix_add_matrix(output4, msk.W0_matrix, jW1, two_k, kss, two_k, kss,order);                                     //TODO Maybe not correct with rho(j) !=0 / =0
                W0_w1_rj = matrix_mul_vector(output5, W0_W1, vj.rj[h].vec_rj, two_k, kss, kss, order);

                for (int s = 0; s < (two_k); ++s) {
                    g2_mul_fix(sk.sk13[h].sk_three[s], t_pre_h, W0_w1_rj[s]);
                }
                h++;
            }

            int gh = 0;
            for (auto it = lsss_row[i].begin(); it != lsss_row[i].end(); ++it) {
                //Computes W_j * rj by matrix-vector multiplication.
                Wr = matrix_mul_vector(output1, msk.W_matrix, vj.rj[gh].vec_rj, two_k, kss, kss, order);
                v_plus_w = vector_add_vector(output1_v_plus_w, vj.vj[gh].vec_j, Wr, two_k, two_k, order);                                                   //TODO Bug (relic stack smashing) here but only occurs when k=1 and N_Attr=50. Every other values work
                //Sets sk_1j by adding all vj vectors with the resulting Wr vectors.
                for (int u = 0; u < (two_k); ++u) {
                    g2_mul_fix(sk.sk13[gh].sk_one[u], t_pre_h, v_plus_w[u]);
                }
                gh++;
            }
        }
    }

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
        gt_null(gt_mul_test); gt_new(gt_mul_test);
        gt_null(gt_st_test); gt_new(gt_st_test);
        fp12_set_dig(gt_st_test, 1);

        //Set M = 1
        gt_null(CT_A.M); gt_new(CT_A.M);
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
    gt_null(exp_val); gt_null(exp_val_extra); gt_null(prod); gt_null(mul_val); gt_null(mul_val_extra);
    gt_new(exp_val); gt_new(exp_val_extra); gt_new(prod); gt_new(mul_val); gt_new(mul_val_extra);

    //Temporary variable supposed to hold intermediate result of the calculations.
    gt_t tmp_res; gt_null(tmp_res); gt_new(tmp_res);

    //Initializes the list of coefficients which should yield a size of N_ATTR * (kss+1)
    int p = 0;
    for (auto it = lsss_row[kss].begin(); it != lsss_row[kss].end(); ++it) {
        bn_null(pack_coef[p]); bn_new(pack_coef[p]);
        p++;
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

        //Using OpenABE recover method to get the (N_ATTR * (kss+1)) coefficients
        lsss_2[0].l_recoverCoefficients(policy, attrList);
        lsss_row[0] = lsss_2[0].l_getRows();

        int j = 0;
        //For all Attributes, set up the two lists used for the pp_map_sim_oatep_k12 operation.
        for (auto it = lsss_row[0].begin(); it != lsss_row[0].end(); ++it) {
            //Copy all the coefficients to the pack_coef list.
            bn_copy(pack_coef[wj], it->second.element().m_ZP);
            //Shitty way of doing this.
            bn_t neg_coef; bn_null(neg_coef); bn_new(neg_coef);
            bn_copy(neg_coef, it->second.element().m_ZP);
            bn_t_negate(neg_coef, order);
            bn_copy(pack_coef_neg[wj], neg_coef);

            //TODO Uncomment these lines below to observe that the coefficients are identical.
            //Set up the two lists used for the pp_map_sim_oatep_k12 operation
            for (int jk = 0; jk < ((2 * two_k) + kss); ++jk) {
                if (jk < two_k) {
                    g1_neg(pair_g1[jk], CT_A.C_1[jk]);
                    g2_copy(pair_g2[jk], sk.sk13[j].sk_one[jk]);
                } else if (jk >= two_k && jk < (2 * two_k)) {
                    g1_neg(pair_g1[jk], CT_A.C_23[j].c_3_vec[jk % two_k]);
                    g2_copy(pair_g2[jk], sk.sk13[j].sk_three[jk % two_k]);
                } else {
                    g1_copy(pair_g1[jk], CT_A.C_23[j].c_2_vec[jk % (2 * two_k)]);
                    g2_copy(pair_g2[jk], sk.sk13[j].sk_two[jk % (2 * two_k)]);
                }
            }

            //TODO not used right now as we only use policies of AND gates.
            for (int ik = 0; ik < (two_k); ++ik) {
                g1_copy(pair_g1_ex[ik], CT_A.C_1[ik]);
                g2_copy(pair_g2_ex[ik], sk.sk4[j].sk_four[ik]);
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
        gt_t final_final_res; gt_null(final_final_res); gt_new(final_final_res);
        gt_mul(final_final_res, tmp_res, CT_A.C_4_one_val);

        //Uncomment for correctness check;
        assert(gt_cmp(final_final_res, CT_A.M) == RLC_EQ);
        std::cout << "[*] PASSED" << std::endl;
    }

    print_results("Results decryption():           ", t, NTESTS);
    printf("]\n");

    //Test if msk and mpk is initialized correctly:
    //print_msk(&msk, N_ATTR, kss);
    //print_mpk(&mpk, N_ATTR, kss);
    //print_sk(&sk, &vj, N_ATTR, kss);
    //print_ct(&CT_A, N_ATTR, kss);
    return 0;
}


