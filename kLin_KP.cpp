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
    std::cout << "Benchmarking KP-ABE from K-Lin\n";
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
    std::cout << keyInputWrong;

    struct master_key_k_lin msk;
    struct public_key_k_lin mpk;

    init_master_key_k_lin(N_ATTR, KSEC, &msk);
    init_public_key_k_lin(N_ATTR, KSEC, &mpk);

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
    g1_t group1;
    g2_t group2;
    g1_null(group1); g1_new(group1);
    g2_null(group2); g2_new(group2);

    g1_rand(group1);
    g2_rand(group2);
    g1_mul_pre(t_pre_g, group1);
    g2_mul_pre(t_pre_h, group2);

    /* Setup */
    for (int jo = 0; jo < NTESTS; jo++) {
        t[jo] = cpucycles();

        bn_t A_tmp[(KSEC + 1) * KSEC];

        //Initializes the v-vector and sets the entries to some random bn_t value modulo the order.
        for (int d = 0; d < ((KSEC + 1) * KSEC); ++d) {
            if (d < (KSEC + 1)){
                bn_rand_mod(msk.v_share[d], order);
            }
            //Initializes the bn_t entries of the A-matrix as just random bn_t value modulo the order.
            //Also initializes the g1 entries of the A-matrix by doing matrix multiplications and sets the A_(i,j) to g1^(AW_(i,j)).
            bn_rand_mod(A_tmp[d], order);
            g1_mul_fix(mpk.a_mat[d], t_pre_g, A_tmp[d]);
        }


        //Define dimensions of matrix A and vector v and calculate the matrix-vector product Av.
        int a_rows = KSEC;
        int a_cols = KSEC + 1;
        int v_rows = KSEC + 1;

        bn_t *Av;
        bn_t output[KSEC];
        Av = matrix_mul_vector(output, A_tmp, msk.v_share, a_rows, a_cols, v_rows, order);
        gt_t map_tmp[KSEC];

        for (int k = 0; k < KSEC; k++) {
            //Initialize the gt entries of the e-mapping matrix doing matrix multiplications exponent-wise.
            pp_map_oatep_k12(map_tmp[k], group1, group2);
            gt_exp(mpk.e_mat[k], map_tmp[k], Av[k]);
        }

        //Initializes the "n" W-matrices (master secret key) by setting every entry in these matrices to some random bn value mod the order.
        //In the K-Lin paper the master secret key consists of w_1,...,w_n and describe w_0 = 0.
        for (int j = 0; j < (N_ATTR + 1); j++) {
            for (int m = 0; m < ((KSEC + 1) * KSEC); m++) {
                if (j != 0) {
                    //Sets entries for Wi where i = 1,...., N_att +1 to random bn_t elements
                    bn_rand_mod(msk.atts[j].w[m], order);
                } else {
                    //Set all entries in W0 to be zero
                    bn_zero(msk.atts[j].w[m]);
                }
            }
            //Here matrix multiplication is being calculated.

            bn_t one_as_bn; bn_null(one_as_bn); bn_new(one_as_bn);
            bn_set_dig(one_as_bn, 1);
            g1_t one_as_g1; g1_null(one_as_g1); g1_new(one_as_g1);
            g1_mul_fix(one_as_g1, t_pre_g, one_as_bn);

            g1_t *AWi;
            g1_t output[KSEC * KSEC];
            AWi = matrixG1_mul_matrixBN(output, mpk.a_mat, msk.atts[j].w, KSEC, (KSEC + 1), (KSEC + 1), KSEC, one_as_g1);

            //Initializes the "n" AW_i (masker public key).
            for (int x = 0; x < (KSEC * KSEC); ++x) {
                g1_copy(mpk.mats[j].w[x], AWi[x]);
            }
        }
    }

    printf("[");
    print_results("Results gen param():           ", t, NTESTS);

    /* Key Generation */
    std::unique_ptr <L_OpenABEFunctionInput> funcInput = nullptr;

    struct secret_key_K_Lin sk;
    struct sk_tmp_vj vj;
    init_secret_key_K_Lin(N_ATTR, &sk);
    init_sk_tmp_vj(N_ATTR, KSEC,&vj);


    //TODO this way of setting up OpenAbe lsss seems to work as expected but maybe I am wrong. Check correctness again.
    //TODO for policies with OR gates this size needs to be changed, according to the K_Lin paper the size would be <=2*N_ATTR
    L_OpenABELSSS *lsss_1 = new L_OpenABELSSS(1);
    lsss_1 = (L_OpenABELSSS *) malloc((N_ATTR * (KSEC + 1)) * sizeof(L_OpenABELSSS));

    L_OpenABELSSSRowMap *lsss_row = new L_OpenABELSSSRowMap;
    lsss_row = (L_OpenABELSSSRowMap *) malloc((N_ATTR * (KSEC + 1)) * sizeof(L_OpenABELSSSRowMap) * sizeof(L_OpenABELSSS));
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

        //For all KSEC+1 secrets in v:
        for (int i = 0; i < (KSEC + 1); ++i) {
            L_ZP s_aux;
            s_aux.isOrderSet = true;
            bn_copy(s_aux.order, order);
            bn_copy(s_aux.m_ZP, msk.v_share[i]);

            //Obtain the j shares for each secret in v. Here we j = N_ATT because the policy consists of only AND gates and so all nodes of the policy tree adds exactly a single share.
            lsss_1[i].l_shareSecret(policy,s_aux);
            lsss_row[i] = lsss_1[i].l_getRows();

            bn_t *Wr;
            bn_t output1[KSEC + 1];
            bn_t output2[KSEC + 1];
            int w_rows = (KSEC + 1); int w_cols = KSEC; int r_rows = KSEC;

            int h = 0;
            for (auto it = lsss_row[i].begin(); it != lsss_row[i].end(); ++it) {
                //Create and set r_j which is a vector of size k of random elements g2 elements, and sets sk_2j = r_j
                for (int k = 0; k < (KSEC); k++) {
                    bn_rand_mod(vj.rj[h].vec_rj[k], order);
                    g2_mul_fix(sk.sk[h].sk_two[k], t_pre_h, vj.rj[h].vec_rj[k]);
                }
                //Sets the vj's to contain the j shares for the (KSEC+1) secrets of v.
                //To clarify each vj is a vector of size (KSEC+1) and there are a total of j vectors.
                bn_copy(vj.vj[h].vec_j[i], it->second.element().m_ZP);
                h++;
            }



            bn_t *v_plus_w;
            bn_t output1_v_plus_w[KSEC + 1];

            int gh = 0;
            for (int kj = 0; kj < N_ATTR; kj++) {
                //Computes W_j * rj by matrix-vector multiplication.
                Wr = matrix_mul_vector(output1, msk.atts[gh+1].w, vj.rj[gh].vec_rj, w_rows, w_cols, r_rows, order);
                v_plus_w = vector_add_vector(output1_v_plus_w, vj.vj[gh].vec_j, Wr, (KSEC + 1), (KSEC + 1), order);


                //Sets sk_1j by adding all vj vectors with the resulting Wr vectors.
                for (int u = 0; u < (KSEC + 1); ++u) {
                    //bn_t_add(output2[u], vj.vj[h].vec_j[u], Wr[u], order);
                    g2_mul_fix(sk.sk[gh].sk_one[u], t_pre_h, v_plus_w[u]);
                }
                gh++;
            }
        }
    }

    print_results("Results keyGen():           ", t, NTESTS);

    /* Encryption */
    //Initialize ciphertext struct
    struct ciphertext_K_Lin CT_A;
    init_ciphertext_K_Lin(N_ATTR, KSEC, &CT_A);
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

        //Sample random vector s of size k and set c_3
        for (int i = 0; i < KSEC; ++i) {
            bn_rand_mod(rnd_s[i], order);
            gt_exp(gt_mul_test, mpk.e_mat[i], rnd_s[i]);
            gt_mul(gt_st_test, gt_st_test, gt_mul_test);
        }

        //Multiply gt value of sTAv with the message M to complete ct_3
        gt_mul(CT_A.C_3_one_val, gt_st_test, CT_A.M);

        //set ct_1
        g1_t *ct_1;
        g1_t output[KSEC + 1];

        //Calculate sT*A using vector-matrix multiplication for a transposed vector.
        ct_1 = vector_trans_mul_matrix_g1(output, rnd_s, mpk.a_mat, KSEC, KSEC+1, KSEC);

        //Finishing ct_1 by doing the exponentiation of g.
        for (int t = 0; t < (KSEC +1); ++t) {
            //g1_mul_fix(CT_A.C_1[t], t_pre_g, ct_1[t]);
            g1_copy(CT_A.C_1[t], ct_1[t]);
        }

        //set ct_2i
        //For all N_ATTR + 1 calculate sTAW_i and the reason for doing it over N_ATTR+1 opposed to N_ATTR is because W_0 = 0 and is done to support that the lsss map can be rho(j) = 0.
        for (int a = 0; a < (N_ATTR + 1); ++a) {
            g1_t *ct2_i;
            g1_t output[KSEC];
            ct2_i = vector_trans_mul_matrix_g1(output, rnd_s, mpk.mats[a].w, KSEC, KSEC, KSEC);

            //Finishing c_2i, by doing the exponentiation of g.
            for (int v = 0; v < KSEC; ++v) {
                //g1_mul_fix(CT_A.C_2[a].c_2_mat[v], t_pre_g, ct2_i[v]);
                g1_copy(CT_A.C_2[a].c_2_mat[v], ct2_i[v]);
            }
        }
    }


    print_results("Results encryption():           ", t, NTESTS);
    //TODO start/complete decryption.

    /* Decryption */
    printf("\n");

    //TODO for policies with OR gates this size needs to be changed, according to the K_Lin paper the size would be <=2*N_ATTR
    //List of all the wj coefficients since j = N_ATTR and because we have (KSEC+1) secrets to be shared from v, the total amount of coefficients becomes N_ATTR * (KSEC+1).
    bn_t pack_coef[N_ATTR * (KSEC+1)];

    //Two lists used for the pp_map_sim_oatep_k12 operation.
    //Since ct2j and sk2j are vectors of size (KSEC) and ct1 and sk1j are vectors of size (KSEC+1) we need to do ((KSEC+1)+KSEC) component wise operations.
    g1_t pair_g1[(KSEC+1) + KSEC];
    g2_t pair_g2[(KSEC+1) + KSEC];

    //Two temporary list that is holding intermediate values of the calculations. The size is related to amount of secret of v and is therefore (KSEC+1).
    gt_t tmp_exp_list[KSEC+1];
    gt_t tmp_mul_list[KSEC+1];

    //Temporary variable supposed to hold intermediate result of the calculations.
    gt_t tmp_res; gt_null(tmp_res); gt_new(tmp_res);

    //Initializes the list of coefficients which should yield a size of N_ATTR * (KSEC+1)
    int p = 0;
    for (int k = 0; k < (KSEC + 1); ++k) {
        for (auto it = lsss_row[k].begin(); it != lsss_row[k].end(); ++it) {
            bn_null(pack_coef[p]);
            bn_new(pack_coef[p]);
            p++;
        }
    }

    for (int go = 0; go < NTESTS; go++) {
        t[go] = cpucycles();

        //TODO should map_sim be reset after every iteration?
        gt_t map_sim; gt_null(map_sim); gt_new(map_sim);

        int wj = 0;

        //Initializes the temporary lists.
        for (int r = 0; r < (KSEC+1); ++r) {
            //printf("\n");
            gt_null(tmp_exp_list[r]);
            gt_new(tmp_exp_list[r]);
            gt_null(tmp_mul_list[r]);
            gt_new(tmp_mul_list[r]);

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
            lsss_1[r].l_recoverCoefficients(policy, attrList);
            lsss_row[r] = lsss_1[r].l_getRows();

            //int sd = r;
            int j = 0;

            //For all Attributes, set up the two lists used for the pp_map_sim_oatep_k12 operation.
            for (auto it = lsss_row[r].begin(); it != lsss_row[r].end(); ++it) {
                //Copy all the coefficients to the pack_coef list.
                bn_copy(pack_coef[wj], it->second.element().m_ZP);

                //TODO Uncomment these lines below to observe that the coefficients are identical.
                //printf("Coefs for secret v_%d: \n", r);
                //bn_print(pack_coef[wj]);

                    //Set up the two lists used for the pp_map_sim_oatep_k12 operation
                    for (int jk = 0; jk < ((KSEC + 1) + KSEC); ++jk) {
                        if (jk < (KSEC + 1)) {
                            g1_neg(pair_g1[jk], CT_A.C_1[jk]);
                            g2_copy(pair_g2[jk], sk.sk[j].sk_one[jk]);
                        } else {
                            g1_copy(pair_g1[jk], CT_A.C_2[j + 1].c_2_mat[(jk+1) % KSEC]);
                            g2_copy(pair_g2[jk], sk.sk[j].sk_two[(jk+1) % KSEC]);
                        }
                    }
                    //TODO maybe exponentiating map_sim with wj (of size k+1) in here and add outside
                    pp_map_sim_oatep_k12(map_sim, pair_g1, pair_g2, ((KSEC + 1) + KSEC));

                    //Not necessary for decryption itself, this calculation is used to test if we can recover the secrets manually using the coefficients. AND WE CAN!
                    bn_t_mul(tt_mul_2, pack_coef[wj], vj.vj[j].vec_j[r], order);
                    bn_t_add(tt_add_2, tt_add_2, tt_mul_2, order);

                    //Here we do map_sim = [-sTAv_j]^(wj) where map_sim = [-sTAv_j] comes from the correctness of the K_Lin paper and wj is the coefficients.
                    gt_exp(tmp_exp_list[r], map_sim,pack_coef[wj]);
                    //Here we basically compute the product of [-sTAv_j]^(wj) and saves the result in tmp_mul_list[r]
                    //TODO seems like the product here is not calculated correctly as we are resetting tmp_mul_list every iteration in the outer loop
                    //TODO despite this, decryption still works?
                    gt_mul(tmp_mul_list[r], tmp_mul_list[r], tmp_exp_list[r]);

                //sd += (KSEC + 1);
                wj++;
                j++;
            }
            //Here we complete the product of [-sTAv_j]^(wj)
            //TODO since we reset tmp_res this product seem to not be fully correct, we are basically just return the last iteration's product.
            //TODO Despite this, decryption still works?
            //TODO Best guess: OpenABE lsss recover coefficients are identical.
            gt_mul(tmp_res, tmp_res, tmp_mul_list[r]);
        }

        //Printouts for correctness.
        gt_t final_final_res;
        gt_null(final_final_res); gt_new(final_final_res);

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

    //Test if msk and mpk is initialized correctly:
    //print_msk(&msk, N_ATTR, KSEC);
    //print_mpk(&mpk, N_ATTR, KSEC);
    //print_sk(&sk, &vj, N_ATTR, KSEC);
    //print_ct(&CT_A, N_ATTR, KSEC);

    return 0;
}
