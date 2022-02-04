#include "relic_type3/relic.h"
#include "relic_type3/relic_test.h"

void ibe();

main() {
    ibe();
    return 0;
}

//This is IBE 1-2-3 ABC Comb my beard
void ibe() {
    core_init();
    g1_t P;
    bn_t s;
    bn_null(s);
    bn_new(s);
    g1_null(P);
    g1_new(P);
    RLC_TRY {
        pc_param_set_any();
        printf("%d", pc_param_level());
        pc_param_print();
        //Generator P
        g1_get_gen(P);
        bn_t whatever;
        bn_null(whatever);
        bn_new(whatever);
        g1_get_ord(whatever);
        type3_bn_rand_mod(s, whatever);
        bn_print(s);

        g1_t p_pub;
        g1_new(p_pub);
        g1_copy(p_pub, P);
        g1_mul_gen(p_pub, s);
        g1_print(p_pub);

        uint8_t id[] = {1};
        uint8_t Q_id[32];
        type3_md_map_sh256(Q_id, id, 1);
        g2_t Q_id_g2;
        g2_null(Q_id_g2);
        g2_new(Q_id_g2);g2_map(Q_id_g2, Q_id, 32);
        g2_t d_id;
        g2_null(d_id);
        g2_new(d_id);
        g2_copy(d_id, Q_id_g2);
        g2_mul(d_id, d_id, s);


        bn_t r;
        bn_null(r);
        bn_new(r);
        bn_rand_mod(r, whatever);

        g1_t C_1;
        g1_null(C_1);
        g1_new(C_1);
        g1_copy(C_1, P);
        g1_mul(C_1, C_1, r);

        //g1_t C_2;
        //g1_null(C_2);
        //g1_new(C_2);

        gt_t g_id;gt_null(g_id);gt_new(g_id);pc_map(g_id, p_pub, Q_id_g2);
        gt_print(g_id);
        gt_exp(g_id, g_id, r);
        gt_print(g_id);
        int buffer_size = gt_size_bin(g_id, 0);
        uint8_t bin[buffer_size];
        gt_write_bin(bin, buffer_size, g_id, 0);
        uint8_t hash[32];
        md_map_sh256(hash, bin, buffer_size);


        uint8_t C_2[32];
        char m[32] = "lol this is pretty easy";
        uint8_t message[32];
        for (int i = 0; i < 32; ++i) {
            message[i] = (uint8_t) m[i];
        }
        //uint8_t message[32] = {1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6,
        //                      7, 8};
        for (int i = 0; i < 32; ++i) {
            C_2[i] = (hash[i] ^ message[i]);
        }

        gt_t d_map;gt_null(d_map);gt_new(d_map);pc_map(d_map, C_1, d_id);
        buffer_size = gt_size_bin(d_map, 0);
        uint8_t dbin[buffer_size];
        gt_write_bin(dbin, buffer_size, d_map, 0);
        uint8_t dhash[32];
        md_map_sh256(dhash, dbin, buffer_size);

        uint8_t decrypted_message[32];

        for (int i = 0; i < 32; ++i) {
            decrypted_message[i] = (dhash[i] ^ C_2[i]);
        }
        char finalm[32];
        for (int i = 0; i < 32; ++i) {
            printf("Index %d of message is %d \n", i, message[i]);
            printf("Index %d of decrypted message is %d \n", i, decrypted_message[i]);
            printf("Index %d of C2 message is %d \n", i, C_2[i]);
            finalm[i] = (char) decrypted_message[i];
        }
        printf("%s", finalm);

    }

    //Test comment
    core_clean();
}
