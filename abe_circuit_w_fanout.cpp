//
// Created by futr on 2/18/22.
//

#include <cstdio>
#include <string>
#include "bench_defs.h"
using namespace std;

int getAttrNumber(string attr) {
    attr.erase(0, std::string("|attr").length() - 1);
    return stoi(attr);
}

uint32_t share(L_OpenABEPolicy policy) {

}



int main(int argc, char **argv){
    int test_attr;
    if (argc < 2){
        test_attr = 5;
        printf("sheesh \n");
    } else {
        test_attr = atoi(argv[1]);
        printf("yeesh \n");
    }

    printf("N_attr = %d", test_attr);
    std::string keyInput = "";
    std::string encInput = "";

    uint32_t N_ATTR = test_attr;

    uint32_t *attr_int_list = NULL;
    attr_int_list = (uint32_t *) malloc(sizeof(uint32_t) * test_attr);

    int d = 1;

    for (int k = 0; k < test_attr; k++) {
        keyInput = keyInput + "attr" + std::to_string(d);
        encInput = encInput + "attr" + std::to_string(d);

        if (k < test_attr - 1) {
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

    printf("yeet\n");


    g1_t g1; g1_null(g1); g1_new(g1);
    g2_t g2; g2_null(g2); g1_new(g2);

    g1_get_gen(g1);
    cout << "g1\n";
    g1_print(g1);
    g2_get_gen(g2);
    cout << "g2\n";
    g2_print(g2);

    gt_t gt; gt_null(gt); gt_new(gt);
    pc_map(gt, g1, g2);
    cout << "gt\n";
    gt_print(gt);



    L_ZP y; bn_null(y.m_ZP); bn_new(y.m_ZP);
    y.isOrderSet = true;
    bn_copy(y.order, order);
    bn_rand_mod(y.m_ZP, order);

    gt_t big_y; gt_null(big_y); gt_new(big_y);
    g1_t big_ts_g1[N_ATTR];
    g2_t big_ts_g2_inv[N_ATTR];

    for (int i = 0; i < N_ATTR; i++) {
        bn_t t_i; bn_null(t_i); bn_new(t_i);
        bn_t t_i_inv; bn_null(t_i_inv); bn_new(t_i_inv);
        bn_rand_mod(t_i, order);
        cout << "t_" << i << "\n";
        bn_print(t_i);
        bn_mod_inv(t_i_inv, t_i, order);
        cout << "1/t_" << i << "\n";
        bn_print(t_i_inv);
        g1_null(big_ts_g1[i]); g1_new(big_ts_g1[i]);
        g2_null(big_ts_g2_inv[i]); g2_new(big_ts_g2_inv[i]);
        g1_mul(big_ts_g1[i], g1, t_i);
        cout << "g1^t_" << i << "\n";
        g1_print(big_ts_g1[i]);
        g2_mul(big_ts_g2_inv[i], g2, t_i_inv);
        cout << "g2^(1/t_" << i << ")\n";
        g2_print(big_ts_g2_inv[i]);
    }
    gt_exp(big_y, gt, y.m_ZP);
    cout << "big_y\n";
    gt_print(big_y);
    //Encryption
    L_OpenABEAttributeList *attrList = nullptr;

    unique_ptr<L_OpenABEFunctionInput> encFuncInput = nullptr;
    encFuncInput = L_createAttributeList(encInput);

    if (encFuncInput == nullptr) {
        printf("Invalid attribute encryption input\n");
        return(EXIT_FAILURE);
    }
    if ((attrList = dynamic_cast<L_OpenABEAttributeList *>(encFuncInput.get())) == nullptr) {
        printf("Error in attribute list\n");
        return(EXIT_FAILURE);
    }

    cout << "Encryption input: " << attrList->l_toString() << "\n";

    g1_t big_es[N_ATTR];
    bn_t s; bn_null(s); bn_new(s)
    bn_rand_mod(s, order);
    cout << "s\n";
    bn_print(s);
    for (int i = 0; i < N_ATTR; i++) {
        g1_null(big_es[i]); g1_new(big_es[i]);
        cout << "Value of g1^t_" << i << "\n";
        g1_print(big_ts_g1[i]);
        g1_copy(big_es[i], big_ts_g1[i]);
        g1_mul(big_es[i], big_es[i], s);
        cout << "Value of big_es[" << i << "]\n";
        g1_print(big_es[i]);
    }
    gt_t e_prime; gt_null(e_prime); gt_new(e_prime);
    //ATM we just encrypt 1
    gt_exp(e_prime, big_y, s);
    //KeyGen
    L_OpenABEPolicy *policy;
    L_OpenABELSSS lsss(1);
    L_OpenABELSSSRowMap lsssRows;

    unique_ptr<L_OpenABEFunctionInput> funcInput = nullptr;


    funcInput = L_createPolicyTree(keyInput);
    if (!funcInput) {
        printf("Create policy error in key generation\n");
        return -1;
    }
    policy = dynamic_cast<L_OpenABEPolicy *>(funcInput.get());
    if (policy == nullptr) {
        printf("Error in input policy\n");
        return -1;
    }


    cout << "Key input: " << policy-> l_toString() << "\n";
    lsss.l_shareSecret(policy, y);
    lsssRows = lsss.l_getRows();
    g2_t big_ds[N_ATTR];
    int attr_int;
    bn_t shares[N_ATTR];
    for (auto it = lsssRows.begin(); it != lsssRows.end(); it++) {
        //My name is Deez
        cout << "Share\n";
        bn_print(it -> second.element().m_ZP);
        string label = it -> second.label();
        attr_int = getAttrNumber(label);
        int attr_index = attr_int-1;
        bn_null(shares[attr_index]); bn_new(shares[attr_index]);
        bn_copy(shares[attr_index], it -> second.element().m_ZP);
        g2_null(big_ds[attr_index); g2_new(big_ds[attr_index]);
        g2_mul(big_ds[attr_index], big_ts_g2_inv[attr_index], shares[attr_index]);
        cout << "D(" << attr_index<< ")\n";
        g2_print(big_ds[attr_index]);
    }



    //Decryption
    gt_t big_vas[N_ATTR];
    lsss.l_recoverCoefficients(policy, attrList);
    lsssRows = lsss.l_getRows();
    g1_t e_reconstruct[N_ATTR];
    bn_t reconCoeffs[N_ATTR];
    for (auto it = lsssRows.begin(); it != lsssRows.end(); it++) {
        gt_null(big_vas[N_ATTR]); gt_new(big_vas[N_ATTR]);
        cout << "recon coefficient\n";
        bn_print(it -> second.element().m_ZP);
        string label = it -> second.label();
        attr_int = getAttrNumber(label);
        int attr_index = attr_int-1;
        bn_null(reconCoeffs[attr_index]); bn_new(reconCoeffs[attr_index]);
        //cout << "recon element " << it -> second.element() << "\n";
        bn_copy(reconCoeffs[attr_index], it -> second.element().m_ZP);
        //cout << "recon element zp " << attr_index << "\n";
        //bn_print(reconCoeffs[attr_index]);
        g1_null(e_reconstruct[attr_index]); g1_new(e_reconstruct[attr_index]);
        g1_mul(e_reconstruct[attr_index], big_es[attr_index], reconCoeffs[attr_index]);
        //cout << "e_reconstruct[" << attr_index << "]\n";
        //g1_print(e_reconstruct[attr_index]);
        //pc_map(big_vas[attr_index], e_reconstruct[attr_index], big_ds[attr_index]);
        //cout << "big_vas[" << attr_index << "]\n";
        //gt_print(big_vas[attr_index]);
    }
    gt_t r; gt_null(r); gt_new(r);
    gt_set_unity(r);
    pc_map_sim(r, e_reconstruct, big_ds, N_ATTR);
    cout << "r\n";
    gt_print(r);

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
    bn_print(y.m_ZP);
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
    cout << "Value of blinding factor before decrypt\n";
    gt_print(e_prime);
    cout << "Value of blinding factor after decrypt\n";
    gt_print(r);








    //printf("Testing CP-ABE context\n");
    //cout << "\tkeyInput: " << keyInput << endl;




    return 0;
}
