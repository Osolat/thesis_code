//
// Created by futr on 2/18/22.
//

#include <cstdio>
#include <string>
#include "bench_defs.h"
using namespace std;

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


int getAttrNumber(string attr) {
    attr.erase(0, std::string("|attr").length() - 1);
    return stoi(attr);
}




void test_abe(uint32_t N_ATTR) {
    printf("abe circuit oe, N_attr = %d", N_ATTR);


    std::string keyInput = "";
    std::string encInput = "";

    uint32_t *attr_int_list = NULL;
    attr_int_list = (uint32_t *) malloc(sizeof(uint32_t) * N_ATTR);

    int d = 1;
    for (int k = 0; k < N_ATTR; k++) {
        keyInput = keyInput + "attr" + std::to_string(d);
        encInput = encInput + "attr" + std::to_string(d);

        if (k < N_ATTR - 1) {
            keyInput = keyInput + " and ";
            encInput = encInput + "|";
        }

        attr_int_list[k] = d;

        d++;

    }
    core_init();

    //Public params

    bn_t order;
    pc_param_set_any();
    pc_param_print();
    pc_get_ord(order);


    g1_t g1; g1_null(g1); g1_new(g1);
    g2_t g2; g2_null(g2); g1_new(g2);

    g1_get_gen(g1);
    //cout << "g1\n";
    //g1_print(g1);
    g2_get_gen(g2);
    //cout << "g2\n";
    //g2_print(g2);

    g1_t t_pre_g1[RLC_EP_TABLE_MAX];
    g2_t t_pre_g2[RLC_EP_TABLE_MAX];
    g1_t t_pre_g1_exp[N_ATTR][RLC_EP_TABLE_MAX];
    for (int i = 0; i < RLC_EP_TABLE; i++) {
        g1_new(t_pre_g1[i]);
        g2_new(t_pre_g2[i]);
        for (int j = 0; j < N_ATTR; j++){
            g1_new(t_pre_g1_exp[j][i])
        }
    }
    bn_t y; bn_null(y); bn_new(y);
    gt_t big_y;
    gt_null(big_y);
    gt_new(big_y);
    g1_t big_ts_g1[N_ATTR];
    bn_t exponents[N_ATTR];
    bn_t inv_exponents[N_ATTR];

    for(int j = 0; j < NTESTS; j++) {
        t[j] = cpucycles();

        g1_mul_pre(t_pre_g1, g1);
        g2_mul_pre(t_pre_g2, g2);


        bn_rand_mod(y, order);


        for (int i = 0; i < N_ATTR; i++) {
            bn_null(exponents[i]); bn_new(exponents[i]);
            bn_rand_mod(exponents[i], order);

            g1_null(big_ts_g1[i]); g1_new(big_ts_g1[i]);
            g1_mul_fix(big_ts_g1[i], t_pre_g1, exponents[i]);
            g1_mul_pre(t_pre_g1_exp[i], big_ts_g1[i]);

            bn_null(inv_exponents[i]); bn_new(inv_exponents[i]);


            //g2_mul_fix(big_ts_g2_inv[i], t_pre_g2, t_i_inv);
            //cout << "t_" << i << "\n";
            //bn_print(t_i);

            //cout << "1/t_" << i << "\n";
            //bn_print(t_i_inv);

            //cout << "g1^t_" << i << "\n";
            //g1_print(big_ts_g1[i]);

            //cout << "g2^(1/t_" << i << ")\n";
            //g2_print(big_ts_g2_inv[i]);
        }
        bn_mod_inv_sim(inv_exponents, exponents, order, N_ATTR);
        pc_map(big_y, g1, g2);
        gt_exp(big_y, big_y, y);
    } printf("["); print_results("Results gen param():           ", t, NTESTS);
    //Encryption


    //cout << "sheesh 1\n";
    bn_t attributes[N_ATTR];
    for (size_t i = 0; i < N_ATTR; i++) {
        bn_null(attributes[i]);
        bn_new(attributes[i]);
        bn_set_dig(attributes[i], i + 1);
    }
    //if ((attrList = dynamic_cast<L_OpenABEAttributeList *>(encFuncInput.get())) == nullptr) {
    //    printf("Error in attribute list\n");
    //    return(EXIT_FAILURE);
    //}



    g1_t big_es[N_ATTR];
    bn_t s; bn_null(s); bn_new(s)
    gt_t e_prime; gt_null(e_prime); gt_new(e_prime);

    //cout << "s\n";
    //bn_print(s);
    for (int j = 0; j < NTESTS; j++) {
        t[j] = cpucycles();
        bn_rand_mod(s, order);
        for (int i = 0; i < N_ATTR; i++) {
            g1_null(big_es[i]); g1_new(big_es[i]);
            g1_mul_fix(big_es[i], t_pre_g1_exp[i], s);
            /*
            cout << "Value of big_es[" << i << "]\n";
            g1_print(big_es[i]);

            cout << "Value of g1^t_" << i << "\n";
            g1_print(big_ts_g1[i]);
            */
        }
        //ATM we just encrypt 1
        gt_exp(e_prime, big_y, s);
    } print_results("Results encryption():           ", t, NTESTS);

    //KeyGen
    //L_OpenABEPolicy *policy;
    //L_OpenABELSSS lsss(1);
    //L_OpenABELSSSRowMap lsssRows;

    //unique_ptr<L_OpenABEFunctionInput> funcInput = nullptr;


    //funcInput = L_createPolicyTree(keyInput);
    //if (!funcInput) {
    //    printf("Create policy error in key generation\n");
    //    return -1;
    //}
    //policy = dynamic_cast<L_OpenABEPolicy *>(funcInput.get());
    //if (policy == nullptr) {
    //    printf("Error in input policy\n");
    //    return -1;
    //}



    struct node tree_root;
    //cout << or_tree_formula(N_ATTR) << "\n";
    //cout << and_tree_formula(N_ATTR) << "\n";
    vector<policy_coefficient> vector;
    tree_from_string(and_tree_formula(N_ATTR), &tree_root);
    //cout << "Key input\n";
    //print_tree(&tree_root);
    //cout << "shuus\n";
    //lsss.l_shareSecret(policy, y);
    //lsssRows = lsss.l_getRows();
    g2_t big_ds[N_ATTR];
    int attr_int;
    bn_t shares[N_ATTR];
    for (int j = 0; j < NTESTS; j++) {
        t[j] = cpucycles();
        vector = std::vector<policy_coefficient>();
        share_secret(&tree_root, y, order, vector, true);

        for (auto it = vector.begin(); it != vector.end(); it++) {

            int attr_index = it -> leaf_index-1;
            bn_t tmp; bn_null(tmp); bn_new(tmp);
            //bn_add(rSecret, rSecret, it -> share);
            //bn_null(shares[attr_index]); bn_new(shares[attr_index]);
            //bn_copy(shares[attr_index], it -> share);

            bn_mul(tmp, inv_exponents[attr_index], it -> share);
            bn_mod(tmp, tmp, order);

            g2_null(big_ds[attr_index]); g2_new(big_ds[attr_index]);

            g2_mul_fix(big_ds[attr_index], t_pre_g2, tmp);
            /*
            cout << "Share\n";
            bn_print(it -> share);

            cout << "D(" << attr_index<< ")\n";
            g2_print(big_ds[attr_index]);

            cout << "leaf index: " << it->leaf_index << "\n";
            */

        }

    } print_results("Results key gen():           ", t, NTESTS);
    //cout << "y\n";
    //bn_print(y);
    //bn_mod(rSecret, rSecret, order);
    //cout << "sum of shares\n";
    //bn_print(rSecret);


    //Decryption
    gt_t big_vas[N_ATTR];
    //lsss.l_recoverCoefficients(policy, attrList);
    //lsssRows = lsss.l_getRows();
    try {
        check_satisfiability(&tree_root, attributes, N_ATTR);
        //std::cout << "Satisfiable with correct attributes" << std::endl;
    } catch (struct TreeUnsatisfiableException *e) {
        std::cout << e->what() << std::endl;
    }
    g1_t e_reconstruct[N_ATTR];
    bn_t reconCoeffs[N_ATTR];

    gt_t r; gt_null(r); gt_new(r);
    for (int j = 0; j < NTESTS; j++) {
        t[j] = cpucycles();
        vector = std::vector<policy_coefficient>();
        vector = recover_coefficients(&tree_root, attributes, N_ATTR);
        for (auto it = vector.begin(); it != vector.end(); it++) {
            g1_null(e_reconstruct[N_ATTR]); g1_new(e_reconstruct[N_ATTR]);

            //string label = it -> second.label();
            //attr_int = getAttrNumber(label);
            int attr_index = it->leaf_index - 1;

            g1_mul(e_reconstruct[attr_index], big_es[attr_index], it -> coeff);

            //cout << "e_reconstruct[" << attr_index << "]\n";
            //g1_print(e_reconstruct[attr_index]);

            //pc_map(big_vas[attr_index], e_reconstruct[attr_index], big_ds[attr_index]);
            //cout << "big_vas[" << attr_index << "]\n";
            //gt_print(big_vas[attr_index]);

            //bn_null(reconCoeffs[attr_index]);bn_new(reconCoeffs[attr_index]);
            //cout << "recon coefficient\n";
            //bn_print(it->coeff);

            //cout << "recon element " << it -> second.element() << "\n";
            //bn_print(reconCoeffs[attr_index]);

            //cout << "leaf index: " << it->leaf_index << "\n";
        }
        pc_map_sim(r, e_reconstruct, big_ds, N_ATTR);
    } print_results("Results decryption():           ", t, NTESTS);
    cout << "]\n";


    /*
    bn_t res; bn_null(res); bn_new(res);
    gt_t gt_res; gt_null(gt_res); gt_new(gt_res);
    gt_t map_res; gt_null(map_res); gt_new(map_res);
    gt_set_unity(gt_res);
    gt_set_unity(map_res);
    for (int i = 0; i < N_ATTR; i++) {
        cout << "Share " << i << "\n";
        bn_print(shares[i]);
        cout << "Recon coefficient " << i << "\n";
        bn_print(reconCoeffs[i]);
        bn_t temp; bn_null(temp); bn_new(temp);
        bn_mul(temp, shares[i], reconCoeffs[i]);
        bn_add(res, res, temp);

        //WHY DOES IT WORK NOW!?!?!?!?!
        //there was an error on line 183 brah
        g1_t a; g1_null(a); g1_new(a);
        g1_copy(a, big_ts_g1[i]);
        g1_mul(a, a, s);
        cout << "is this big_e[" << i << "]?: " << (g1_cmp(a, big_es[i]) == RLC_EQ) << "\n";
        g1_mul(a, a, reconCoeffs[i]);
        cout << "is this e_reconstructs?[" << i << "]?: " << (g1_cmp(a, e_reconstruct[i]) == RLC_EQ) << "\n";
        g2_t b; g2_null(b); g2_new(b);
        g2_copy(b, big_ts_g2_inv[i]);
        g2_mul(b, b, shares[i]);
        cout << "is this big_ds[" << i << "]?: " << (g2_cmp(b, big_ds[i]) == RLC_EQ) << "\n";
        if (g2_cmp(b, big_ds[i]) != RLC_EQ) {
            g2_print(b);
            cout << "*************************************\n";
            g2_print(big_ds[i]);
        }
        gt_t stuff1; gt_null(stuff1); gt_new(stuff1);
        pc_map(stuff1, a, b);
        //gt_exp(stuff1, stuff1, s);
        cout << "e(g1^t_" << i << "share[" << i << "], g2^t_" << i << "coeff[" << i << "])s\n";
        gt_print(stuff1);
        gt_t stuff2; gt_null(stuff2); gt_new(stuff2);
        gt_exp(stuff2, gt, temp);
        gt_exp(stuff2, stuff2, s);
        cout << "gt^(share[" << i << "]*coeff[" << i << "])s\n";
        gt_mul(map_res, map_res, stuff1);
        gt_mul(gt_res, gt_res, stuff2);
    }
    bn_mod(res, res, order);
    cout << "Secret before recon\n";
    bn_print(y);
    cout << "Secret after recon\n";
    bn_print(res);
    cout << "pow(gt,res*s)\n";
    gt_t x2; gt_null(x2); gt_new(x2);
    gt_exp(x2, gt, res);
    gt_exp(x2, x2, s);
    gt_print(x2);
    cout << "gt_res\n";
    gt_print(gt_res);
    cout << "map_res\n";
    gt_print(map_res);
    //pp_map_sim_oatep_k12(blindingFactor, big_es, d_reconstruct, N_ATTR);
     */
    //
    cout << "[*] Correctness check r = E': " << (gt_cmp(r, e_prime) == RLC_EQ) << endl;

    cout << "Value of blinding factor before decrypt\n";
    gt_print(e_prime);
    cout << "Value of blinding factor after decrypt\n";
    gt_print(r);
}

int main(int argc, char **argv) {
    t[0] = cpucycles();
    int test_attr;


    uint32_t N_ATTR = test_attr;

    //uint32_t *attr_int_list = NULL;
    //attr_int_list = (uint32_t *) malloc(sizeof(uint32_t) * test_attr);
    test_abe(2);
    //test_abe(8);
    //test_abe(16);
    return 0;
}