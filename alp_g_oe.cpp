#include <cstdio>
#include <string>
#include <cmath>
#include "bench_defs.h"
#include "alp_common.cpp"
using namespace std;

void test_g(int N_ATTR) {
    int bound = N_ATTR+1;
    printf("alp naive, N_attr = %d", N_ATTR);


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

    //cout << "Public params\n";

    bn_t order;
    pc_param_set_any();
    pc_param_print();
    pc_get_ord(order);

    
    cout << "[";
    struct alp_pp_naive_oe pp;
    bn_t alpha; bn_null(alpha); bn_new(alpha);
    bn_rand_mod(alpha, order);
    setup_naive_oe(&pp, alpha, order, bound);
    //print_public_params(pp, bound);

    //cout << "Key Gen\n";
    struct alp_sk_oe sk;
    struct node tree_root;
    tree_from_string(and_tree_formula(N_ATTR), &tree_root);
    lsss_vector = std::vector<policy_coefficient>();
    init_secret_key_alp_oe(bound, &sk);
    bn_t shares[N_ATTR]; bn_t r[N_ATTR];
    for (size_t j = 0; j < NTESTS; j++) {
        t[j] = cpucycles();
        keygen_naive_oe(pp, &sk, &tree_root, alpha);
    } print_results("Results KeyGen():          ", t, NTESTS);
    //print_secret_key_oe(sk, bound); 


    bn_t attributes[N_ATTR];
    bn_t p_Coeffs[bound];
    for (size_t i = 0; i < N_ATTR; i++) {
        bn_null(attributes[i]); 
        bn_new(attributes[i]); 
        bn_set_dig(attributes[i], i+1);
    }
    
    struct alp_ciphertext_oe C;
    for (size_t j = 0; j < NTESTS; j++){
        t[j] = cpucycles();
        coeff_array(p_Coeffs, attributes, bound);
        encrypt_naive_oe(pp, p_Coeffs, &C);
    } print_results("Results Encrypt():          ", t, NTESTS);

    try {
        check_satisfiability(&tree_root, attributes, N_ATTR);
        //std::cout << "Satisfiable with correct attributes" << std::endl;
    } catch (struct TreeUnsatisfiableException *e) {
        std::cout << e->what() << std::endl;
    }

    for (size_t j = 0; j < NTESTS; j++) {
        t[j] = cpucycles();
        decrypt_naive_oe(pp, sk, C, attributes, tree_root);
    } print_results("Results Decrypt():         ", t, NTESTS);
    cout << "]\n"; 
    free_tree(&tree_root);
}


int main (int argc, char **argv) {
    test_g(2);
    test_g(4);
    test_g(8);
    test_g(16);
    test_g(32);
    test_g(64);
    test_g(128);
    test_g(256);
    test_g(512);
    test_g(1024);
    return 0;
}