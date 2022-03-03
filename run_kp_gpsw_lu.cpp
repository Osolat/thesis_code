//
// Created by benjamin on 2/11/22.
//

#include <algorithm>
#include <cstdio>
#include <string>

#include "bench_defs.h"

unsigned char *uint32_to_u_char_array(const uint32_t n) {
    unsigned char *a;

    a = (unsigned char *)malloc(4 * sizeof(unsigned char));
    a[0] = (n >> 24) & 0xff;
    a[1] = (n >> 16) & 0xff;
    a[2] = (n >> 8) & 0xff;
    a[3] = n & 0xff;

    return a;
}
long long cpucycles(void) {
    unsigned long long result;
    asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
                 : "=a"(result)::"%rdx");
    return result;
}

static int cmp_llu(const void *a, const void *b) {
    if (*(unsigned long long *)a < *(unsigned long long *)b)
        return -1;
    if (*(unsigned long long *)a > *(unsigned long long *)b)
        return 1;
    return 0;
}

static unsigned long long median(unsigned long long *l, size_t llen) {
    qsort(l, llen, sizeof(unsigned long long), cmp_llu);

    if (llen % 2)
        return l[llen / 2];
    else
        return (l[llen / 2 - 1] + l[llen / 2]) / 2;
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
    if (argc == 1) {
        printf("Need to give argument\n");
        return 0;
    }

    int test_attr = atoi(argv[1]);

    std::string keyInput = "";
    std::string encInput = "";

    uint32_t N_ATTR = test_attr;

    uint32_t *attr_int_list = NULL;
    attr_int_list = (uint32_t *)malloc(sizeof(uint32_t) * test_attr);

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

    struct master_key_kp_gpsw_lu msk;
    struct public_key_kp_gpsw_lu mpk;

    init_master_key_kp_gpsw_lu(N_ATTR, &msk);
    init_public_key_kp_gpsw_lu(N_ATTR, &mpk);

    core_init();

    bn_t order;
    pc_param_set_any();
    pc_param_print();
    pc_get_ord(order);

    /* Setup */

    std::unique_ptr<L_OpenABEFunctionInput> keyFuncInput = nullptr;
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

    /*y random in ZP*/
    bn_rand_mod(msk.y, order);

    /*g1 = g^y*/
    g1_mul_gen(mpk.g1, msk.y);

    /*g2 random in G2*/
    g2_rand(mpk.g2);

    /*KeyGeneration*/
    struct secret_key_kp_gpsw_lu sk;
    init_secret_key_kp_gpsw_lu(test_attr, &sk);

    std::unique_ptr<L_OpenABEFunctionInput> funcInput = nullptr;
    L_OpenABELSSS lsss(1);
    L_OpenABELSSSRowMap lsssRows;
    L_OpenABEPolicy *policy;

    funcInput = L_createPolicyTree(encInput);
    if (!funcInput) {
        printf("Create policy error in KeyGeneration\n");
        return -1;
    }

    policy = dynamic_cast<L_OpenABEPolicy *>(funcInput.get());

    if (policy == nullptr) {
        printf("Error in input policy\n");
        return -1;
    }

    /*Secret sharing of y, according to policy tree*/
    L_ZP s_aux;
    s_aux.isOrderSet = true;
    bn_copy(s_aux.order, order);
    bn_copy(s_aux.m_ZP, msk.y);
    lsss.l_shareSecret(policy, s_aux);
    lsssRows = lsss.l_getRows();

    g2_t temp;
    g2_new(temp);
    g2_null(temp);

    bn_t r;
    bn_new(r);
    bn_null(r);
    /*Accessing q_leaf(0) <= second.element().m_ZP*/
    int i = 0;
    uint8_t hash[32];
    for (auto it = lsssRows.begin(); it != lsssRows.end(); ++it) {
        /*Dx = g_2^(q_x(0)) T(i)^(r_x)*/
        // TODO: Precomputation tables
        g2_mul(sk.D_values[i], mpk.g2, it->second.element().m_ZP);
        /*  g2_mul(sk.D_values[i], mpk.g2, it->second.element().m_ZP);
         std::string attribute = it->second.label();
         std::copy(attribute.begin(), attribute.end(), std::begin(hash)); */
        // Has form attr-----
        // TODO: Confirm this is the correct way of doing T(i) aka hashing
        std::string attr = it->second.label();
        attr.erase(0, std::string("|attr").length() - 1);
        uint32_t attr_int_1 = stoi(attr);

        unsigned char *input_data = uint32_to_u_char_array(attr_int_1);
        g2_map(temp, input_data, 4);
        // md_map_sh256(hash, hash, 32);
        // g2_map(temp, hash, 32);
        std::cout << attr_int_1 << std::endl;
        g2_print(temp);
        // Temp is the hash of attribute string into G2 at this point
        bn_rand_mod(r, order);
        g2_mul(temp, temp, r);
        // g2^qx * T(i)^r
        g2_add(sk.D_values[i], sk.D_values[i], temp);
        g1_mul_gen(sk.R_values[i], r);
        i++;
    }

    /* Encryption */
    gt_t message;
    gt_new(message);
    gt_null(message);
    gt_rand(message);

    bn_t s;
    bn_new(s);
    bn_null(s);
    bn_rand_mod(s, order);
    struct ciphertext_kp_gpsw_lu E;
    init_ciphertext_kp_gpsw_lu(test_attr, &E);
    /*E = (gamma, me(g1,g2)^s, g^s, E_i = T(i)^s)*/
    /*g^s*/
    g1_mul_gen(E.E_prime_prime, s);

    /*m * e(g1,g2)^s*/
    gt_t e;
    gt_new(e);
    gt_null(e);
    pp_map_oatep_k12(e, mpk.g1, mpk.g2);
    gt_exp(e, e, s);
    gt_mul(E.E_prime, message, e);
    i = 0;
    for (auto it = lsssRows.begin(); it != lsssRows.end(); ++it) {
        /* std::string attribute = it->second.label();
        // Has form attr-----
        std::copy(attribute.begin(), attribute.end(), std::begin(hash));
        // Has form attr-----
        md_map_sh256(hash, hash, 32);
        g2_map(temp, hash, 32); */
        // TODO: Confirm this is the correct way of doing T(i) aka hashing
        std::string attr = it->second.label();
        attr.erase(0, std::string("|attr").length() - 1);
        uint32_t attr_int_1 = stoi(attr);
        unsigned char *input_data = uint32_to_u_char_array(attr_int_1);
        g2_map(temp, input_data, 4);
        std::cout << attr_int_1 << std::endl;
        g2_print(temp);
        // Temp is the hash of attribute string into G2 at this point
        g2_mul(E.E_values[i], temp, s);
        i++;
    }

    /*Decryption*/
    lsss.l_recoverCoefficients(policy, attrList);
    lsssRows = lsss.l_getRows();

    gt_t F_root;
    gt_new(F_root);
    gt_null(F_root);
    // TODO: Is this legal? fp12_set_dig
    fp12_set_dig(F_root, 1);
    gt_t mapping;
    gt_new(mapping);
    gt_null(mapping);

    gt_t mapping2;
    gt_new(mapping2);
    gt_null(mapping2);

    i = 0;
    for (auto it = lsssRows.begin(); it != lsssRows.end(); ++it) {
        // TODO: Fix so it works with special attribute policies
        // TODO: Could do pp_map_sim here instead. Sim is better
        pp_map_oatep_k12(mapping, E.E_prime_prime, sk.D_values[i]);
        pp_map_oatep_k12(mapping2, sk.R_values[i], E.E_values[i]);
        gt_inv(mapping2, mapping2);
        gt_mul(mapping, mapping, mapping2);
        gt_exp(mapping, mapping, it->second.element().m_ZP);
        gt_mul(F_root, F_root, mapping);
        i++;
    }

    gt_t result;
    gt_new(result);
    gt_null(result);

    gt_inv(F_root, F_root);
    gt_mul(result, F_root, E.E_prime);

    printf("------------------ \n");
    gt_print(result);
    printf("----------------------\n");
    printf("------------------ \n");
    gt_print(message);
    printf("----------------------\n");
    printf("Result of comparison between Message and F_root: %d\n", gt_cmp(message, result) == RLC_EQ);
    return 0;
}
