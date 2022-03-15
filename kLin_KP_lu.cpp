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

    init_master_key_k_lin_lu(N_ATTR, KSEC, &msk);
    init_public_key_k_lin_lu(N_ATTR, KSEC, &mpk);

    core_init();

    bn_t order;
    pc_param_set_any();
    pc_param_print();
    g1_get_ord(order);

    //DEBUG UTIL:
    //test_matrix_mul_vector(KSEC, order);
    //test_vector_trans_mul_matrix(KSEC, order);
    //test_matrix_mul_matrix(KSEC, order);
    //test_vector_dot_product(KSEC, order);
    //test_matrix_mul_scalar(KSEC, order);
    //test_matrix_add_matrix(KSEC, order);
    //test_vector_add_vector(KSEC, order);

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

    int two_k = ((2*KSEC) + 1);

    /* Setup */
    for (int jo = 0; jo < NTESTS; jo++) {
        t[jo] = cpucycles();

        bn_t A1_tmp[two_k * KSEC];


        for (int d = 0; d < (two_k * KSEC); ++d) {
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
        bn_t output[KSEC];
        g1_t outputAW[KSEC * KSEC];
        g1_t outputAW0[KSEC * KSEC];
        g1_t outputAW1[KSEC * KSEC];
        gt_t map_tmp[KSEC];

        bn_t one_as_bn;
        bn_null(one_as_bn);
        bn_new(one_as_bn);
        bn_set_dig(one_as_bn, 1);
        g1_t one_as_g1;
        g1_null(one_as_g1);
        g1_new(one_as_g1);
        g1_mul_fix(one_as_g1, t_pre_g, one_as_bn);

        Av = matrix_mul_vector(output, A1_tmp, msk.v_secret, KSEC, two_k, two_k, order);
        AW = matrixG1_mul_matrixBN(outputAW, mpk.A1_mat, msk.W_matrix, KSEC, two_k, two_k, KSEC, one_as_g1);
        AW0 = matrixG1_mul_matrixBN(outputAW0, mpk.A1_mat, msk.W0_matrix, KSEC, two_k, two_k, KSEC,one_as_g1);
        AW1 = matrixG1_mul_matrixBN(outputAW1, mpk.A1_mat, msk.W1_matrix, KSEC, two_k, two_k, KSEC, one_as_g1);

        for (int k = 0; k < (KSEC * KSEC); k++) {
            if (k < KSEC) {
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
    init_sk_tmp_vectors_lu(N_ATTR, KSEC, &vj);


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
            bn_t output3[two_k * KSEC];
            bn_t output4[two_k * KSEC];
            bn_t output5[two_k];
            bn_t output1_v_plus_w[two_k];

            int h = 0;
            for (auto it = lsss_row[i].begin(); it != lsss_row[i].end(); ++it) {
                //Create and set r_j which is a vector of size k of random elements g2 elements, and sets sk_2j = r_j
                for (int k = 0; k < (KSEC); k++) {
                    bn_rand_mod(vj.rj[h].vec_rj[k], order);
                    g2_mul_fix(sk.sk13[h].sk_two[k], t_pre_h, vj.rj[h].vec_rj[k]);
                }
                //Sets the vj's to contain the j shares for the (KSEC+1) secrets of v.
                //To clarify each vj is a vector of size (KSEC+1) and there are a total of j vectors.
                bn_copy(vj.vj[h].vec_j[i], it->second.element().m_ZP);
                g2_mul_fix(sk.sk4[h].sk_four[i], t_pre_h, vj.vj[h].vec_j[i]);                                                                                                         //Correct if rho(j)=0 for all j, however should be empty when all AND gates in policy tree

                jW1 = matrix_mul_scalar(output3, msk.W1_matrix, (h), two_k, KSEC, order);                                                               //h+1 so that we don't multiply with 0
                W0_W1 = matrix_add_matrix(output4, msk.W0_matrix, jW1, two_k, KSEC, two_k, KSEC,order);                                     //TODO Maybe not correct with rho(j) !=0 / =0
                W0_w1_rj = matrix_mul_vector(output5, W0_W1, vj.rj[h].vec_rj, two_k, KSEC, KSEC, order);

                for (int s = 0; s < (two_k); ++s) {
                    g2_mul_fix(sk.sk13[h].sk_three[s], t_pre_h, W0_w1_rj[s]);
                }
                h++;
            }

            int gh = 0;
            for (auto it = lsss_row[i].begin(); it != lsss_row[i].end(); ++it) {
                //Computes W_j * rj by matrix-vector multiplication.
                Wr = matrix_mul_vector(output1, msk.W_matrix, vj.rj[gh].vec_rj, two_k, KSEC, KSEC, order);
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
    init_ciphertext_K_Lin_lu(N_ATTR, KSEC, &CT_A);
    init_tmp_si_lu(N_ATTR, KSEC, &si);
    bn_t rnd_s[KSEC];

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
        for (int i = 0; i < KSEC; ++i) {
            bn_rand_mod(rnd_s[i], order);
            gt_exp(gt_mul_test, mpk.e_mat[i], rnd_s[i]);
            gt_mul(gt_st_test, gt_st_test, gt_mul_test);
        }

        //Multiply gt value of sTAv with the message M to complete ct_4
        gt_mul(CT_A.C_4_one_val, gt_st_test, CT_A.M);


        //set ct_1
        g1_t *ct_1_g1;
        g1_t output_g1[two_k];
        ct_1_g1 = vector_trans_mul_matrix_g1(output_g1, rnd_s, mpk.A1_mat, KSEC, two_k, KSEC);

        for (int t = 0; t < (two_k); ++t) {
            g1_copy(CT_A.C_1[t], ct_1_g1[t]);
        }

        for (int z = 0; z < N_ATTR; ++z) {                                                                                                                                              //i=N_ATTR because all attribute is needed to decrypt due the fact it's all AND gates
            for (int x = 0; x < KSEC; ++x) {
                bn_rand_mod(si.si[z].si_vec[x], order);
            }

            //set ct_3
            g1_t *ct_3;
            g1_t output_ct_3[two_k];

            //Calculate sT*A using vector-matrix multiplication for a transposed vector.
            ct_3 = vector_trans_mul_matrix_g1(output_ct_3, si.si[z].si_vec, mpk.A1_mat, KSEC, two_k, KSEC);

            //Finishing ct_1 by doing the exponentiation of g.
            for (int ti = 0; ti < (two_k); ++ti) {
                g1_copy(CT_A.C_23[z].c_3_vec[ti], ct_3[ti]);
            }

            g1_t *sTAW;
            g1_t output[KSEC];
            sTAW = vector_trans_mul_matrix_g1(output, rnd_s, mpk.AW_mat, KSEC, KSEC, KSEC);

            g1_t *i_A1_W1; g1_t *W0_i_A1_W1; g1_t *res_sT_W0; g1_t *res_Add_val;
            g1_t output_i_A1_W1[KSEC * KSEC];
            g1_t output_W0_i_A1_W1[KSEC * KSEC];
            g1_t output_res_sT_W0[KSEC];
            g1_t output_res_Add_val[KSEC];

            i_A1_W1 = matrix_mul_scalar_g1(output_i_A1_W1, mpk.AW1_mat, (z), KSEC, KSEC);
            W0_i_A1_W1 = matrix_add_matrix_g1(output_W0_i_A1_W1, mpk.AW0_mat, i_A1_W1, KSEC, KSEC, KSEC,KSEC);
            res_sT_W0 = vector_trans_mul_matrix_g1(output_res_sT_W0, si.si[z].si_vec, W0_i_A1_W1, KSEC, KSEC, KSEC);
            res_Add_val = vector_add_vector_g1(output_res_Add_val, sTAW, res_sT_W0, KSEC, KSEC);

            for (int p = 0; p < KSEC; ++p) {
                g1_copy(CT_A.C_23[z].c_2_vec[p], res_Add_val[p]);
            }
        }
    }


    print_results("Results encryption():           ", t, NTESTS);

    //TODO start/complete decryption.
    printf("\n");

    /* Decryption */
    bn_t pack_coef[N_ATTR * two_k];
    bn_t pack_coef_neg[N_ATTR * two_k];

    g1_t pair_g1[(2 * two_k) + KSEC];
    g2_t pair_g2[(2 * two_k) + KSEC];
    g1_t pair_g1_ex[two_k];
    g2_t pair_g2_ex[two_k];

    //Two temporary list that is holding intermediate values of the calculations. The size is related to amount of secret of v and is therefore (KSEC+1).
    gt_t tmp_exp_list[two_k];
    gt_t tmp_mul_list[two_k];
    gt_t tmp_prod_list[two_k];
    gt_t tmp_exp_list_extra[two_k];
    gt_t tmp_mul_list_extra[two_k];

    //Temporary variable supposed to hold intermediate result of the calculations.
    gt_t tmp_res; gt_null(tmp_res); gt_new(tmp_res);

    //Initializes the list of coefficients which should yield a size of N_ATTR * (KSEC+1)
    int p = 0;
    for (int k = 0; k < (two_k); ++k) {
        for (auto it = lsss_row[k].begin(); it != lsss_row[k].end(); ++it) {
            bn_null(pack_coef[p]); bn_new(pack_coef[p]);
            p++;
        }
    }

    for (int go = 0; go < NTESTS; go++) {
        t[go] = cpucycles();

        //TODO should map_sim be reset after every iteration?
        gt_t map_sim; gt_null(map_sim); gt_new(map_sim);
        gt_t map_sim_2; gt_null(map_sim_2); gt_new(map_sim_2);

        int wj = 0;
        for (int r = 0; r < (two_k); ++r) {
            gt_null(tmp_exp_list[r]); gt_new(tmp_exp_list[r]);
            gt_null(tmp_mul_list[r]); gt_new(tmp_mul_list[r]);
            gt_null(tmp_exp_list_extra[r]); gt_new(tmp_exp_list_extra[r]);
            gt_null(tmp_mul_list_extra[r]); gt_new(tmp_mul_list_extra[r]);
            gt_null(tmp_prod_list[r]); gt_new(tmp_prod_list[r]);

            //Sets tmp_mul_list[r] to one so that the multiplication starts out correct.
            //TODO if not here decryption fails, however seems wrong to reset to one for each iteration of the outer loop.
            fp12_set_dig(tmp_mul_list[r], 1);

            bn_t tt_mul_2; bn_null(tt_mul_2); bn_new(tt_mul_2);
            bn_t tt_add_2; bn_null(tt_add_2); bn_new(tt_add_2);

            //Sets tmp_res to one so that the final multiplications starts out correct.
            //TODO if not here decryption fails, however seems wrong to reset to one for each iteration since we then dont get the full product and dont use all coefficients.
            fp12_set_dig(tmp_res, 1);

            //Using OpenABE recover method to get the (N_ATTR * (KSEC+1)) coefficients
            //TODO this way of using OpenABE lsss/recoverCoefficients seems to work but some strange things is observed here.
            //TODO Observation: The first j=N_ATTR recovered coefficients are identical to the next j=N_ATTR recovered coefficients, and again they are identical the the next set of coefficients.
            //TODO Seems to be a OpenABE recoverCoefficients related issue as no matter the original secret and different shares the coefficients always seems identical.
            lsss_2[r].l_recoverCoefficients(policy, attrList);
            lsss_row[r] = lsss_2[r].l_getRows();

            //int sd = r;
            int j = 0;
            //For all Attributes, set up the two lists used for the pp_map_sim_oatep_k12 operation.
            for (auto it = lsss_row[r].begin(); it != lsss_row[r].end(); ++it) {
                //Copy all the coefficients to the pack_coef list.
                bn_copy(pack_coef[wj], it->second.element().m_ZP);
                //Shitty way of doing this.
                bn_t neg_coef; bn_null(neg_coef); bn_new(neg_coef);
                bn_copy(neg_coef, it->second.element().m_ZP);
                bn_t_negate(neg_coef, order);
                bn_copy(pack_coef_neg[wj], neg_coef);

                //TODO Uncomment these lines below to observe that the coefficients are identical.
                //Set up the two lists used for the pp_map_sim_oatep_k12 operation
                for (int jk = 0; jk < ((2 * two_k) + KSEC); ++jk) {
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
                for (int ik = 0; ik < (two_k); ++ik) {
                    g1_copy(pair_g1_ex[ik], CT_A.C_1[ik]);
                    g2_copy(pair_g2_ex[ik], sk.sk4[j].sk_four[ik]);
                }

                //TODO maybe exponentiating map_sim with wj (of size k+1) in here and add outside
                //printf("Do this amount of times \n");
                pp_map_sim_oatep_k12(map_sim, pair_g1, pair_g2, ((2 * two_k) + KSEC));
                pp_map_sim_oatep_k12(map_sim_2, pair_g1_ex, pair_g2_ex, two_k);

                //Not necessary for decryption itself, this calculation is used to test if we can recover the secrets manually using the coefficients. AND WE CAN!
                bn_t_mul(tt_mul_2, pack_coef[wj], vj.vj[j].vec_j[r], order);
                bn_t_add(tt_add_2, tt_add_2, tt_mul_2, order);

                //Here we do map_sim = [-sTAv_j]^(wj) where map_sim = [-sTAv_j] comes from the correctness of the K_Lin paper and wj is the coefficients.
                gt_exp(tmp_exp_list[r], map_sim, pack_coef[wj]);

                //TODO use for rho(j)=0
                gt_exp(tmp_exp_list_extra[r], map_sim_2, pack_coef_neg[wj]);

                gt_mul(tmp_prod_list[r], tmp_exp_list[r], tmp_exp_list_extra[r]);

                gt_mul(tmp_mul_list[r], tmp_mul_list[r], tmp_exp_list[r]);

                //TODO multiply cases where X_rho(j)=1 and rho(j)=0
                //sd += (KSEC + 1);
                wj++;
                j++;
            }
            //Here we complete the product of [-sTAv_j]^(wj)
            //TODO since we reset tmp_res this product seem to not be fully correct, we are basically just return the last iteration's product.
            //TODO Despite this, decryption still works?
            gt_mul(tmp_res, tmp_res, tmp_mul_list[r]);
        }

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
    //print_msk(&msk, N_ATTR, KSEC);
    //print_mpk(&mpk, N_ATTR, KSEC);
    //print_sk(&sk, &vj, N_ATTR, KSEC);
    //print_ct(&CT_A, N_ATTR, KSEC);
    return 0;
}


