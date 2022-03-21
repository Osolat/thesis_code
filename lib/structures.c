/*
 * This file is part of the ABE Squared (https://github.com/abecryptools/abe_squared).
 * Copyright (c) 2022 Antonio de la Piedra, Marloes Venema and Greg Alpár
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "structures.h"

#include <stdlib.h>

#include "g1_arith.h"
#include "g2_arith.h"
#include "gt_arith.h"
#include "zp_arith.h"

// rw13 oe
int init_ciphertext_rw13_ok(const uint32_t n_attr, struct ciphertext_rw13_ok *c) {
    c->N_ATTR = n_attr;

    c->C_1 = (struct c_attribute_rw13_ok_1 *)malloc(n_attr * sizeof(struct c_attribute_rw13_ok_1));
    c->C_2 = (struct c_attribute_rw13_ok_2 *)malloc(n_attr * sizeof(struct c_attribute_rw13_ok_2));
    c->C_3 = (struct c_attribute_rw13_ok_3 *)malloc(n_attr * sizeof(struct c_attribute_rw13_ok_3));

    if (c->C_1 == NULL || c->C_2 == NULL || c->C_3 == NULL) {
        return EXIT_FAILURE;
    } else {
        g2_null(&c->C_PRIMA);
        g2_new(&c->C_PRIMA);

        for (uint32_t i = 0; i < n_attr; i++) {
            g2_null(&c->C_1[i].c_attr);
            g2_new(&c->C_1[i].c_attr);

            g2_null(&c->C_2[i].c_attr);
            g2_new(&c->C_2[i].c_attr);

            g2_null(&c->C_3[i].c_attr);
            g2_new(&c->C_3[i].c_attr);
        }

        return EXIT_SUCCESS;
    }
}

int init_ciphertext_rw13_oe(const uint32_t n_attr, struct ciphertext_rw13_oe *c) {
    c->N_ATTR = n_attr;

    c->C_1 = (struct c_attribute_rw13_oe_1 *)malloc(n_attr * sizeof(struct c_attribute_rw13_oe_1));
    c->C_2 = (struct c_attribute_rw13_oe_2 *)malloc(n_attr * sizeof(struct c_attribute_rw13_oe_2));
    c->C_3 = (struct c_attribute_rw13_oe_3 *)malloc(n_attr * sizeof(struct c_attribute_rw13_oe_3));

    if (c->C_1 == NULL || c->C_2 == NULL || c->C_3 == NULL) {
        return EXIT_FAILURE;
    } else {
        g1_null(&c->C_PRIMA);
        g1_new(&c->C_PRIMA);

        for (uint32_t i = 0; i < n_attr; i++) {
            g1_null(&c->C_1[i].c_attr);
            g1_new(&c->C_1[i].c_attr);

            g1_null(&c->C_2[i].c_attr);
            g1_new(&c->C_2[i].c_attr);

            g1_null(&c->C_3[i].c_attr);
            g1_new(&c->C_3[i].c_attr);
        }

        return EXIT_SUCCESS;
    }
}

void get_k_attr_rw13_ok(const uint32_t attr, const struct secret_key_rw13_ok s, g1_t k_attr_1, g1_t k_attr_2) {
    for (uint32_t i = 0; i < s.N_ATTR; i++) {
        if (s.attributes_K_1[i].attr == attr) {
            g1_copy(k_attr_1, s.attributes_K_1[i].k_attr);
            g1_copy(k_attr_2, s.attributes_K_2[i].k_attr);
        }
    }
}

void get_k_attr_rw13_oe(const uint32_t attr, const struct secret_key_rw13_oe s, g2_t k_attr_1, g2_t k_attr_2) {
    for (uint32_t i = 0; i < s.N_ATTR; i++) {
        if (s.attributes_K_1[i].attr == attr) {
            g2_copy(k_attr_1, s.attributes_K_1[i].k_attr);
            g2_copy(k_attr_2, s.attributes_K_2[i].k_attr);
        }
    }
}

int init_secret_key_rw13_ok(const uint32_t n_attr, struct secret_key_rw13_ok *s) {
    s->N_ATTR = n_attr;

    s->attributes_K_1 = (struct k_attribute_K_1_ok *)malloc(n_attr * sizeof(struct k_attribute_K_1_ok));
    s->attributes_K_2 = (struct k_attribute_K_2_ok *)malloc(n_attr * sizeof(struct k_attribute_K_2_ok));

    if (s->attributes_K_1 == NULL || s->attributes_K_2 == NULL) {
        return EXIT_FAILURE;
    } else {
        g1_null(&s->K);
        g1_new(&s->K)
            g1_null(&s->K_PRIMA);
        g1_new(&s->K_PRIMA);

        for (uint32_t i = 0; i < n_attr; i++) {
            g1_null(&s->attributes_K_1[i].k_attr);
            g1_null(&s->attributes_K_2[i].k_attr);
            g1_new(&s->attributes_K_1[i].k_attr);
            g1_new(&s->attributes_K_2[i].k_attr);
        }

        return EXIT_SUCCESS;
    }
}

int init_secret_key_rw13_oe(const uint32_t n_attr, struct secret_key_rw13_oe *s) {
    s->N_ATTR = n_attr;

    s->attributes_K_1 = (struct k_attribute_K_1_oe *)malloc(n_attr * sizeof(struct k_attribute_K_1_oe));
    s->attributes_K_2 = (struct k_attribute_K_2_oe *)malloc(n_attr * sizeof(struct k_attribute_K_2_oe));

    if (s->attributes_K_1 == NULL || s->attributes_K_2 == NULL) {
        return EXIT_FAILURE;
    } else {
        g2_null(&s->K);
        g2_new(&s->K)
            g2_null(&s->K_PRIMA);
        g2_new(&s->K_PRIMA);

        for (uint32_t i = 0; i < n_attr; i++) {
            g2_null(&s->attributes_K_1[i].k_attr);
            g2_null(&s->attributes_K_2[i].k_attr);
            g2_new(&s->attributes_K_1[i].k_attr);
            g2_new(&s->attributes_K_2[i].k_attr);
        }

        return EXIT_SUCCESS;
    }
}

int init_public_key_rw13_ok(const uint32_t n_attr, struct public_key_rw13_ok *p) {
    p->N_ATTR = n_attr;
    g1_null(&p->g);
    g1_new(&p->g);
    g2_null(&p->h);
    g2_new(&p->h);
    g2_null(&p->g_b);
    g2_new(&p->g_b);
    g2_null(&p->g_b_0);
    g2_new(&p->g_b_0);
    g2_null(&p->g_b_1);
    g2_new(&p->g_b_1);
    g2_null(&p->g_b_prima);
    g2_new(&p->g_b_prima);

    return EXIT_SUCCESS;
}

int init_public_key_rw13_oe(const uint32_t n_attr, struct public_key_rw13_oe *p) {
    p->N_ATTR = n_attr;
    g1_null(&p->g);
    g1_new(&p->g);
    g2_null(&p->h);
    g2_new(&p->h);
    g1_null(&p->g_b);
    g1_new(&p->g_b);
    g1_null(&p->g_b_0);
    g1_new(&p->g_b_0);
    g1_null(&p->g_b_1);
    g1_new(&p->g_b_1);
    g1_null(&p->g_b_prima);
    g1_new(&p->g_b_prima);

    return EXIT_SUCCESS;
}

int init_public_key_rw13_cp(const uint32_t n_attr, struct public_key_rw13_cp *p) {
    p->N_ATTR = n_attr;
    p->attributes = (struct attribute_blind_wat11_i_ok *)malloc(n_attr * sizeof(struct attribute_blind_wat11_i_ok));

    if (p->attributes == NULL) {
        return EXIT_FAILURE;
    } else {
        g1_null(&p->g2);
        g1_new(&p->g2);
        g2_null(&p->g);
        g2_new(&p->g);
        g1_null(&p->u);
        g1_new(&p->u);
        g1_null(&p->v);
        g1_new(&p->v);
        g1_null(&p->h);
        g1_new(&p->h);
        g1_null(&p->w);
        g1_new(&p->w);
        return EXIT_SUCCESS;
    }
}

int init_ciphertext_rw13_cp(const uint32_t n_attr, struct ciphertext_rw13_cp *c) {
    c->N_ATTR = n_attr;

    c->C_1 = (struct c_attribute_rw13_1 *)malloc(n_attr * sizeof(struct c_attribute_rw13_1));
    c->C_2 = (struct c_attribute_rw13_2 *)malloc(n_attr * sizeof(struct c_attribute_rw13_2));
    c->C_3 = (struct c_attribute_rw13_3 *)malloc(n_attr * sizeof(struct c_attribute_rw13_3));

    if (c->C_1 == NULL || c->C_2 == NULL || c->C_3 == NULL) {
        return EXIT_FAILURE;
    } else {
        g2_null(&c->C_0);
        g2_new(&c->C_0);

        for (uint32_t i = 0; i < n_attr; i++) {
            g1_null(&c->C_1[i].c_attr);
            g1_new(&c->C_1[i].c_attr);

            g1_null(&c->C_2[i].c_attr);
            g1_new(&c->C_2[i].c_attr);

            g2_null(&c->C_3[i].c_attr);
            g2_new(&c->C_3[i].c_attr);
        }

        return EXIT_SUCCESS;
    }
}

int init_secret_key_rw13_cp(const uint32_t n_attr, struct secret_key_rw13_cp *s) {
    s->N_ATTR = n_attr;

    s->attributes_K_2 = (struct k_attribute_K_2 *)malloc(n_attr * sizeof(struct k_attribute_K_2));
    s->attributes_K_3 = (struct k_attribute_K_3 *)malloc(n_attr * sizeof(struct k_attribute_K_3));

    if (s->attributes_K_3 == NULL || s->attributes_K_2 == NULL) {
        return EXIT_FAILURE;
    } else {
        g1_null(&s->K_0);
        g1_new(&s->K_0)
            g2_null(&s->K_1);
        g2_new(&s->K_1)

            for (uint32_t i = 0; i < n_attr; i++) {
            g2_null(&s->attributes_K_2[i].k_attr);
            g1_null(&s->attributes_K_3[i].k_attr);
            g2_new(&s->attributes_K_2[i].k_attr);
            g1_new(&s->attributes_K_3[i].k_attr);
        }

        return EXIT_SUCCESS;
    }
}

void get_k_attr_rw13_cp(const uint32_t attr, const struct secret_key_rw13_cp s, g2_t k_attr_2, g1_t k_attr_3) {
    for (uint32_t i = 0; i < s.N_ATTR; i++) {
        if (s.attributes_K_2[i].attr == attr) {
            g2_copy(k_attr_2, s.attributes_K_2[i].k_attr);
            g1_copy(k_attr_3, s.attributes_K_3[i].k_attr);
        }
    }
}

// wat11 ok

void get_c_attr_c_1_wat11_i_ok(const uint32_t attr, const struct ciphertext_wat11_i_ok c, g2_t c_1) {
    for (uint32_t i = 0; i < c.N_ATTR; i++) {
        if (c.C_1[i].attr == attr)
            g2_copy(c_1, c.C_1[i].c_attr);
    }
}

void get_c_attr_c_2_wat11_i_ok(const uint32_t attr, const struct ciphertext_wat11_i_ok c, g2_t c_2) {
    for (uint32_t i = 0; i < c.N_ATTR; i++) {
        if (c.C_2[i].attr == attr)
            g2_copy(c_2, c.C_2[i].c_attr);
    }
}

int init_ciphertext_wat11_i_ok(const uint32_t n_attr, struct ciphertext_wat11_i_ok *c) {
    c->N_ATTR = n_attr;

    c->C_1 = (struct c_attribute_wat11_i_ok *)malloc(n_attr * sizeof(struct c_attribute_wat11_i_ok));
    c->C_2 = (struct c_attribute_wat11_i_ok *)malloc(n_attr * sizeof(struct c_attribute_wat11_i_ok));

    if (c->C_1 == NULL || c->C_2 == NULL) {
        return EXIT_FAILURE;
    } else {
        g2_null(&c->C_PRIMA);
        g2_new(&c->C_PRIMA);

        for (uint32_t i = 0; i < n_attr; i++) {
            g2_null(&c->C_1[i].c_attr);
            g2_new(&c->C_1[i].c_attr);

            g2_null(&c->C_2[i].c_attr);
            g2_new(&c->C_2[i].c_attr);
        }

        return EXIT_SUCCESS;
    }
}

int init_ciphertext_wat11_iv_ok(const uint32_t n_attr, struct ciphertext_wat11_iv_ok *c) {
    c->N_ATTR = n_attr;

    c->C_1 = (struct c_attribute_wat11_iv_ok_1 *)malloc(n_attr * sizeof(struct c_attribute_wat11_iv_ok_1));
    c->C_2 = (struct c_attribute_wat11_iv_ok_2 *)malloc(n_attr * sizeof(struct c_attribute_wat11_iv_ok_2));

    if (c->C_1 == NULL || c->C_2 == NULL) {
        return EXIT_FAILURE;
    } else {
        g2_null(&c->C_PRIMA);
        g2_new(&c->C_PRIMA);

        for (uint32_t i = 0; i < n_attr; i++) {
            g1_null(&c->C_1[i].c_attr);
            g1_new(&c->C_1[i].c_attr);

            g2_null(&c->C_2[i].c_attr);
            g2_new(&c->C_2[i].c_attr);
        }

        return EXIT_SUCCESS;
    }
}

void get_k_attr_wat11_i_ok(const uint32_t attr, const struct secret_key_wat11_i_ok s, g1_t k_attr) {
    for (uint32_t i = 0; i < s.N_ATTR; i++) {
        if (s.attributes[i].attr == attr)
            g1_copy(k_attr, s.attributes[i].k_attr);
    }
}

void get_k_attr_wat11_iv_ok(const uint32_t attr, const struct secret_key_wat11_iv_ok s, g1_t k_attr) {
    for (uint32_t i = 0; i < s.N_ATTR; i++) {
        if (s.attributes[i].attr == attr)
            g1_copy(k_attr, s.attributes[i].k_attr);
    }
}

int init_secret_key_wat11_i_ok(const uint32_t n_attr, struct secret_key_wat11_i_ok *s) {
    s->N_ATTR = n_attr;
    s->attributes = (struct k_attribute_wat11_i_ok *)malloc(n_attr * sizeof(struct k_attribute_wat11_i_ok));

    if (s->attributes == NULL) {
        return EXIT_FAILURE;
    } else {
        g1_null(&s->K_PRIMA);
        g1_new(&s->K_PRIMA)
            g1_null(&s->K);
        g1_new(&s->K)

            for (uint32_t i = 0; i < n_attr; i++) {
            g1_null(&s->attributes[i].k_attr);
            g1_new(&s->attributes[i].k_attr);
        }

        return EXIT_SUCCESS;
    }
}

int init_secret_key_wat11_iv_ok(const uint32_t n_attr, struct secret_key_wat11_iv_ok *s) {
    s->N_ATTR = n_attr;
    s->attributes = (struct k_attribute_wat11_i_ok *)malloc(n_attr * sizeof(struct k_attribute_wat11_i_ok));

    if (s->attributes == NULL) {
        return EXIT_FAILURE;
    } else {
        g1_null(&s->K_PRIMA);
        g2_new(&s->K_PRIMA)
            g2_null(&s->K);
        g1_new(&s->K)

            for (uint32_t i = 0; i < n_attr; i++) {
            g1_null(&s->attributes[i].k_attr);
            g1_new(&s->attributes[i].k_attr);
        }

        return EXIT_SUCCESS;
    }
}

int init_secret_key_wat11_iv_oe(const uint32_t n_attr, struct secret_key_wat11_iv_oe *s) {
    s->N_ATTR = n_attr;
    s->attributes = (struct k_attribute_wat11_iv_oe *)malloc(n_attr * sizeof(struct k_attribute_wat11_iv_oe));

    if (s->attributes == NULL) {
        return EXIT_FAILURE;
    } else {
        g2_null(&s->K_PRIMA);
        g2_new(&s->K_PRIMA)
            g2_null(&s->K);
        g2_new(&s->K)

            for (uint32_t i = 0; i < n_attr; i++) {
            g1_new(&s->attributes[i].k_attr);
            g1_new(&s->attributes[i].k_attr);
        }

        return EXIT_SUCCESS;
    }
}

void get_b_attr_blind_wat11_i_ok(const uint32_t attr, const struct public_key_wat11_i_ok p, g2_t g2_o) {
    for (uint32_t i = 0; i < p.N_ATTR; i++) {
        if (p.attributes[i].attr == attr)
            g2_copy(g2_o, p.attributes[i].h_b_attr);
    }
}

int init_public_key_wat11_i_ok(const uint32_t n_attr, struct public_key_wat11_i_ok *p) {
    p->N_ATTR = n_attr;
    p->attributes = (struct attribute_blind_wat11_i_ok *)malloc(n_attr * sizeof(struct attribute_blind_wat11_i_ok));

    if (p->attributes == NULL) {
        return EXIT_FAILURE;
    } else {
        g1_null(&p->g);
        g1_new(&p->g);
        g2_null(&p->B);
        g2_new(&p->B);
        g2_null(&p->h);
        g2_new(&p->h);

        for (uint32_t i = 0; i < n_attr; i++) {
            g2_null(&p->attributes[i].h_b_attr);
            g2_new(&p->attributes[i].h_b_attr);
        }
        return EXIT_SUCCESS;
    }
}

// wat11 ok

int init_public_key(const uint32_t n_attr, struct public_key *p) {
    p->N_ATTR = n_attr;
    p->attributes = (struct attribute_blind *)malloc(n_attr * sizeof(struct attribute_blind));

    if (p->attributes == NULL) {
        return EXIT_FAILURE;
    } else {
        g1_ar_init(&p->g);
        g1_ar_init(&p->B);
        g2_ar_init(&p->h);

        for (uint32_t i = 0; i < n_attr; i++)
            g1_ar_init(&p->attributes[i].g_b_attr);

        return EXIT_SUCCESS;
    }
}

int init_master_key(const uint32_t n_attr, struct master_key *m) {
    m->N_ATTR = n_attr;
    m->attributes = (struct attribute *)malloc(n_attr * sizeof(struct attribute));

    if (m->attributes == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_ciphertext_wat11_i_oe(const uint32_t n_attr, struct ciphertext_wat11_i_oe *c) {
    c->N_ATTR = n_attr;

    c->C_1 = (struct c_attribute_wat11_i_oe *)malloc(n_attr * sizeof(struct c_attribute_wat11_i_oe));
    c->C_2 = (struct c_attribute_wat11_i_oe *)malloc(n_attr * sizeof(struct c_attribute_wat11_i_oe));

    if (c->C_1 == NULL || c->C_2 == NULL) {
        return EXIT_FAILURE;
    } else {
        g1_null(&c->C_PRIMA);
        g1_new(&c->C_PRIMA);

        for (uint32_t i = 0; i < n_attr; i++) {
            g1_null(&c->C_1[i].c_attr);
            g1_new(&c->C_1[i].c_attr);

            g1_null(&c->C_2[i].c_attr);
            g1_new(&c->C_2[i].c_attr);
        }

        return EXIT_SUCCESS;
    }
}

int init_ciphertext_wat11_iv_oe(const uint32_t n_attr, struct ciphertext_wat11_iv_oe *c) {
    c->N_ATTR = n_attr;

    c->C_1 = (struct c_attribute_wat11_iv_oe_1 *)malloc(n_attr * sizeof(struct c_attribute_wat11_iv_oe_1));
    c->C_2 = (struct c_attribute_wat11_iv_oe_2 *)malloc(n_attr * sizeof(struct c_attribute_wat11_iv_oe_2));

    if (c->C_1 == NULL || c->C_2 == NULL) {
        return EXIT_FAILURE;
    } else {
        g1_null(&c->C_PRIMA);
        g1_new(&c->C_PRIMA);

        for (uint32_t i = 0; i < n_attr; i++) {
            g1_null(&c->C_1[i].c_attr);
            g1_new(&c->C_1[i].c_attr);

            g2_null(&c->C_2[i].c_attr);
            g2_new(&c->C_2[i].c_attr);
        }

        return EXIT_SUCCESS;
    }
}

int init_public_key_wat11_i_oe(const uint32_t n_attr, struct public_key_wat11_i_oe *p) {
    p->N_ATTR = n_attr;
    p->attributes = (struct attribute_blind_wat11_i_oe *)malloc(n_attr * sizeof(struct attribute_blind_wat11_i_oe));

    if (p->attributes == NULL) {
        return EXIT_FAILURE;
    } else {
        g1_null(&p->g);
        g1_new(&p->g);
        g1_null(&p->B);
        g1_new(&p->B);
        g2_null(&p->h);
        g2_new(&p->h);

        for (uint32_t i = 0; i < n_attr; i++) {
            g1_null(&p->attributes[i].g_b_attr);
            g1_new(&p->attributes[i].g_b_attr);
        }

        return EXIT_SUCCESS;
    }
}

int init_public_key_ok(const uint32_t n_attr, struct public_key_ok *p) {
    p->N_ATTR = n_attr;
    p->attributes = (struct attribute_blind_ok *)malloc(n_attr * sizeof(struct attribute_blind_ok));

    if (p->attributes == NULL) {
        return EXIT_FAILURE;
    } else {
        g1_ar_init(&p->g);
        g2_ar_init(&p->B);
        g2_ar_init(&p->h);

        for (uint32_t i = 0; i < n_attr; i++)
            g2_ar_init(&p->attributes[i].h_b_attr);

        return EXIT_SUCCESS;
    }
}

int init_public_key_lu_ok(const uint32_t n_attr, struct public_key_lu_ok *p) {
    p->N_ATTR = n_attr;
    p->attributes = (struct attribute_blind *)malloc(n_attr * sizeof(struct attribute_blind));

    if (p->attributes == NULL) {
        return EXIT_FAILURE;
    } else {
        g1_ar_init(&p->g);
        g1_ar_init(&p->B);
        g2_ar_init(&p->h);

        for (uint32_t i = 0; i < n_attr; i++)
            g1_ar_init(&p->attributes[i].g_b_attr);

        return EXIT_SUCCESS;
    }
}

int init_public_key_cp(const uint32_t n_attr, struct public_key_cp *p) {
    p->N_ATTR = n_attr;
    p->attributes = (struct attribute_blind *)malloc(n_attr * sizeof(struct attribute_blind));

    if (p->attributes == NULL) {
        return EXIT_FAILURE;
    } else {
        g1_ar_init(&p->g);
        g1_ar_init(&p->B);
        g2_ar_init(&p->h);

        for (uint32_t i = 0; i < n_attr; i++)
            g1_ar_init(&p->attributes[i].g_b_attr);

        return EXIT_SUCCESS;
    }
}

int init_secret_key(const uint32_t n_attr, struct secret_key *s) {
    s->N_ATTR = n_attr;
    s->attributes = (struct k_attribute *)malloc(n_attr * sizeof(struct k_attribute));

    if (s->attributes == NULL) {
        return EXIT_FAILURE;
    } else {
        g2_ar_init(&s->K_PRIMA);
        g2_ar_init(&s->K);

        for (uint32_t i = 0; i < n_attr; i++)
            g2_ar_init(&s->attributes[i].k_attr);

        return EXIT_SUCCESS;
    }
}

int init_secret_key_lu_oe(const uint32_t n_attr, struct secret_key_lu_oe *s) {
    s->N_ATTR = n_attr;
    s->attributes = (struct k_attribute_lu_oe *)malloc(n_attr * sizeof(struct k_attribute_lu_oe));

    if (s->attributes == NULL) {
        return EXIT_FAILURE;
    } else {
        g2_ar_init(&s->K_PRIMA);
        g2_ar_init(&s->K);

        for (uint32_t i = 0; i < n_attr; i++)
            g1_ar_init(&s->attributes[i].k_attr);

        return EXIT_SUCCESS;
    }
}

int init_secret_key_ok(const uint32_t n_attr, struct secret_key_ok *s) {
    s->N_ATTR = n_attr;
    s->attributes = (struct k_attribute_ok *)malloc(n_attr * sizeof(struct k_attribute_ok));

    if (s->attributes == NULL) {
        return EXIT_FAILURE;
    } else {
        g1_ar_init(&s->K_PRIMA);
        g1_ar_init(&s->K);

        for (uint32_t i = 0; i < n_attr; i++)
            g1_ar_init(&s->attributes[i].k_attr);

        return EXIT_SUCCESS;
    }
}

int init_secret_key_lu_ok(const uint32_t n_attr, struct secret_key_lu_ok *s) {
    s->N_ATTR = n_attr;
    s->attributes = (struct k_attribute_lu_ok *)malloc(n_attr * sizeof(struct k_attribute_lu_ok));

    if (s->attributes == NULL) {
        return EXIT_FAILURE;
    } else {
        g2_ar_init(&s->K_PRIMA);
        g1_ar_init(&s->K);

        for (uint32_t i = 0; i < n_attr; i++)
            g1_ar_init(&s->attributes[i].k_attr);

        return EXIT_SUCCESS;
    }
}

int init_secret_key_od(const uint32_t n_attr, struct secret_key_od *s) {
    s->N_ATTR = n_attr;
    s->attributes = (struct k_attribute_od *)malloc(n_attr * sizeof(struct k_attribute_od));

    if (s->attributes == NULL) {
        return EXIT_FAILURE;
    } else {
        g2_ar_init(&s->K_PRIMA);
        g2_ar_init(&s->K);

        for (uint32_t i = 0; i < n_attr; i++)
            g1_ar_init(&s->attributes[i].k_attr);

        return EXIT_SUCCESS;
    }
}

int init_secret_key_cp(const uint32_t n_attr, struct secret_key_cp *s) {
    s->N_ATTR = n_attr;
    s->attributes = (struct k_attribute_cp *)malloc(n_attr * sizeof(struct k_attribute_cp));

    if (s->attributes == NULL) {
        return EXIT_FAILURE;
    } else {
        g2_ar_init(&s->K_PRIMA);
        g1_ar_init(&s->K);

        for (uint32_t i = 0; i < n_attr; i++)
            g1_ar_init(&s->attributes[i].k_attr);

        return EXIT_SUCCESS;
    }
}

int init_ciphertext(const uint32_t n_attr, struct ciphertext *c) {
    c->N_ATTR = n_attr;

    c->C_1 = (struct c_attribute *)malloc(n_attr * sizeof(struct c_attribute));
    c->C_2 = (struct c_attribute *)malloc(n_attr * sizeof(struct c_attribute));

    if (c->C_1 == NULL || c->C_2 == NULL) {
        return EXIT_FAILURE;
    } else {
        g1_ar_init(&c->C_PRIMA);

        for (uint32_t i = 0; i < n_attr; i++) {
            g1_ar_init(&c->C_1[i].c_attr);
            g1_ar_init(&c->C_2[i].c_attr);
        }

        return EXIT_SUCCESS;
    }
}

int init_ciphertext_ok(const uint32_t n_attr, struct ciphertext_ok *c) {
    c->N_ATTR = n_attr;

    c->C_1 = (struct c_attribute_ok *)malloc(n_attr * sizeof(struct c_attribute_ok));
    c->C_2 = (struct c_attribute_ok *)malloc(n_attr * sizeof(struct c_attribute_ok));

    if (c->C_1 == NULL || c->C_2 == NULL) {
        return EXIT_FAILURE;
    } else {
        g2_ar_init(&c->C_PRIMA);

        for (uint32_t i = 0; i < n_attr; i++) {
            g2_ar_init(&c->C_1[i].c_attr);
            g2_ar_init(&c->C_2[i].c_attr);
        }

        return EXIT_SUCCESS;
    }
}

int init_ciphertext_od(const uint32_t n_attr, struct ciphertext_od *c) {
    c->N_ATTR = n_attr;

    c->C_1 = (struct c_attribute *)malloc(n_attr * sizeof(struct c_attribute));
    c->C_2 = (struct c_attribute_od *)malloc(n_attr * sizeof(struct c_attribute_od));

    if (c->C_1 == NULL || c->C_2 == NULL) {
        return EXIT_FAILURE;
    } else {
        g1_ar_init(&c->C_PRIMA);

        for (uint32_t i = 0; i < n_attr; i++) {
            g1_ar_init(&c->C_1[i].c_attr);
            g2_ar_init(&c->C_2[i].c_attr);
        }

        return EXIT_SUCCESS;
    }
}

int init_ciphertext_cp(const uint32_t n_attr, struct ciphertext_cp *c) {
    c->N_ATTR = n_attr;

    c->C_1 = (struct c_attribute *)malloc(n_attr * sizeof(struct c_attribute));
    c->C_2 = (struct c_attribute_cp *)malloc(n_attr * sizeof(struct c_attribute_cp));

    if (c->C_1 == NULL || c->C_2 == NULL) {
        return EXIT_FAILURE;
    } else {
        g2_ar_init(&c->C_PRIMA);

        for (uint32_t i = 0; i < n_attr; i++) {
            g1_ar_init(&c->C_1[i].c_attr);
            g2_ar_init(&c->C_2[i].c_attr);
        }

        return EXIT_SUCCESS;
    }
}

int init_ciphertext_lu_oe(const uint32_t n_attr, struct ciphertext_lu_oe *c) {
    c->N_ATTR = n_attr;

    c->C_1 = (struct c_attribute_lu_oe_1 *)malloc(n_attr * sizeof(struct c_attribute_lu_oe_1));
    c->C_2 = (struct c_attribute_lu_oe_2 *)malloc(n_attr * sizeof(struct c_attribute_lu_oe_2));

    if (c->C_1 == NULL || c->C_2 == NULL) {
        return EXIT_FAILURE;
    } else {
        g1_ar_init(&c->C_PRIMA);

        for (uint32_t i = 0; i < n_attr; i++) {
            g1_ar_init(&c->C_1[i].c_attr);
            g2_ar_init(&c->C_2[i].c_attr);
        }

        return EXIT_SUCCESS;
    }
}

int init_ciphertext_lu_ok(const uint32_t n_attr, struct ciphertext_lu_ok *c) {
    c->N_ATTR = n_attr;

    c->C_1 = (struct c_attribute_lu_ok_1 *)malloc(n_attr * sizeof(struct c_attribute_lu_ok_1));
    c->C_2 = (struct c_attribute_lu_ok_2 *)malloc(n_attr * sizeof(struct c_attribute_lu_ok_2));

    if (c->C_1 == NULL || c->C_2 == NULL) {
        return EXIT_FAILURE;
    } else {
        g2_ar_init(&c->C_PRIMA);

        for (uint32_t i = 0; i < n_attr; i++) {
            g1_ar_init(&c->C_1[i].c_attr);
            g2_ar_init(&c->C_2[i].c_attr);
        }

        return EXIT_SUCCESS;
    }
}

void get_b_attr(const uint32_t attr, const struct master_key m, bn_t bn_o) {
    for (uint32_t i = 0; i < m.N_ATTR; i++) {
        if (m.attributes[i].attr == attr)
            bn_t_copy(bn_o, m.attributes[i].b_attr);
    }
}

void set_b_attr(const uint32_t attr, bn_t b_attr, struct master_key *m) {
    for (uint32_t i = 0; i < m->N_ATTR; i++) {
        if (m->attributes[i].attr == attr)
            bn_t_copy(m->attributes[i].b_attr, b_attr);
    }
}

void get_b_attr_blind(const uint32_t attr, const struct public_key p, g1_t g1_o) {
    for (uint32_t i = 0; i < p.N_ATTR; i++) {
        if (p.attributes[i].attr == attr)
            g1_copy(g1_o, p.attributes[i].g_b_attr);
    }
}

void get_b_attr_blind_ok(const uint32_t attr, const struct public_key_ok p, g2_t g2_o) {
    for (uint32_t i = 0; i < p.N_ATTR; i++) {
        if (p.attributes[i].attr == attr)
            g2_copy(g2_o, p.attributes[i].h_b_attr);
    }
}

void get_b_attr_blind_cp(const uint32_t attr, const struct public_key_cp p, g1_t g1_o) {
    for (uint32_t i = 0; i < p.N_ATTR; i++) {
        if (p.attributes[i].attr == attr)
            g1_copy(g1_o, p.attributes[i].g_b_attr);
    }
}

void set_b_attr_blind(const uint32_t attr, g1_t g_b_attr, struct public_key *p) {
    for (uint32_t i = 0; i < p->N_ATTR; i++) {
        if (p->attributes[i].attr == attr)
            g1_copy(p->attributes[i].g_b_attr, g_b_attr);
    }
}

void set_b_attr_blind_ok(const uint32_t attr, g2_t h_b_attr, struct public_key_ok *p) {
    for (uint32_t i = 0; i < p->N_ATTR; i++) {
        if (p->attributes[i].attr == attr)
            g2_copy(p->attributes[i].h_b_attr, h_b_attr);
    }
}

void set_b_attr_blind_cp(const uint32_t attr, g1_t g_b_attr, struct public_key_cp *p) {
    for (uint32_t i = 0; i < p->N_ATTR; i++) {
        if (p->attributes[i].attr == attr)
            g1_copy(p->attributes[i].g_b_attr, g_b_attr);
    }
}

void get_k_attr(const uint32_t attr, const struct secret_key s, g2_t k_attr) {
    for (uint32_t i = 0; i < s.N_ATTR; i++) {
        if (s.attributes[i].attr == attr)
            g2_copy(k_attr, s.attributes[i].k_attr);
    }
}

void get_k_attr_lu_oe(const uint32_t attr, const struct secret_key_lu_oe s, g1_t k_attr) {
    for (uint32_t i = 0; i < s.N_ATTR; i++) {
        if (s.attributes[i].attr == attr)
            g1_copy(k_attr, s.attributes[i].k_attr);
    }
}

void get_k_attr_ok(const uint32_t attr, const struct secret_key_ok s, g1_t k_attr) {
    for (uint32_t i = 0; i < s.N_ATTR; i++) {
        if (s.attributes[i].attr == attr)
            g1_copy(k_attr, s.attributes[i].k_attr);
    }
}

void get_k_attr_lu_ok(const uint32_t attr, const struct secret_key_lu_ok s, g1_t k_attr) {
    for (uint32_t i = 0; i < s.N_ATTR; i++) {
        if (s.attributes[i].attr == attr)
            g1_copy(k_attr, s.attributes[i].k_attr);
    }
}

void get_k_attr_od(const uint32_t attr, const struct secret_key_od s, g1_t k_attr) {
    for (uint32_t i = 0; i < s.N_ATTR; i++) {
        if (s.attributes[i].attr == attr)
            g1_copy(k_attr, s.attributes[i].k_attr);
    }
}

void get_k_attr_cp(const uint32_t attr, const struct secret_key_cp s, g1_t k_attr) {
    for (uint32_t i = 0; i < s.N_ATTR; i++) {
        if (s.attributes[i].attr == attr)
            g1_copy(k_attr, s.attributes[i].k_attr);
    }
}

void get_k_attr_wat11_iv_oe(const uint32_t attr, const struct secret_key_wat11_iv_oe s, g1_t k_attr) {
    for (uint32_t i = 0; i < s.N_ATTR; i++) {
        if (s.attributes[i].attr == attr)
            g1_copy(k_attr, s.attributes[i].k_attr);
    }
}

void set_k_attr(const uint32_t attr, g2_t k_attr, struct secret_key *s) {
    for (uint32_t i = 0; i < s->N_ATTR; i++) {
        if (s->attributes[i].attr == attr)
            g2_copy(s->attributes[i].k_attr, k_attr);
    }
}

void set_k_attr_lu_oe(const uint32_t attr, g1_t k_attr, struct secret_key_lu_oe *s) {
    for (uint32_t i = 0; i < s->N_ATTR; i++) {
        if (s->attributes[i].attr == attr)
            g1_copy(s->attributes[i].k_attr, k_attr);
    }
}

void set_k_attr_ok(const uint32_t attr, g1_t k_attr, struct secret_key_ok *s) {
    for (uint32_t i = 0; i < s->N_ATTR; i++) {
        if (s->attributes[i].attr == attr)
            g1_copy(s->attributes[i].k_attr, k_attr);
    }
}

void set_k_attr_lu_ok(const uint32_t attr, g1_t k_attr, struct secret_key_lu_ok *s) {
    for (uint32_t i = 0; i < s->N_ATTR; i++) {
        if (s->attributes[i].attr == attr)
            g1_copy(s->attributes[i].k_attr, k_attr);
    }
}

void set_k_attr_od(const uint32_t attr, g1_t k_attr, struct secret_key_od *s) {
    for (uint32_t i = 0; i < s->N_ATTR; i++) {
        if (s->attributes[i].attr == attr)
            g1_copy(s->attributes[i].k_attr, k_attr);
    }
}

void set_k_attr_cp(const uint32_t attr, g1_t k_attr, struct secret_key_cp *s) {
    for (uint32_t i = 0; i < s->N_ATTR; i++) {
        if (s->attributes[i].attr == attr)
            g1_copy(s->attributes[i].k_attr, k_attr);
    }
}

void get_c_attr_c_1(const uint32_t attr, const struct ciphertext c, g1_t c_1) {
    for (uint32_t i = 0; i < c.N_ATTR; i++) {
        if (c.C_1[i].attr == attr)
            g1_copy(c_1, c.C_1[i].c_attr);
    }
}

void get_c_attr_c_2(const uint32_t attr, const struct ciphertext c, g1_t c_2) {
    for (uint32_t i = 0; i < c.N_ATTR; i++) {
        if (c.C_2[i].attr == attr)
            g1_copy(c_2, c.C_2[i].c_attr);
    }
}

void get_c_attr_c_1_lu_oe(const uint32_t attr, const struct ciphertext_lu_oe c, g1_t c_1) {
    for (uint32_t i = 0; i < c.N_ATTR; i++) {
        if (c.C_1[i].attr == attr)
            g1_copy(c_1, c.C_1[i].c_attr);
    }
}

void get_c_attr_c_1_lu_ok(const uint32_t attr, const struct ciphertext_lu_ok c, g1_t c_1) {
    for (uint32_t i = 0; i < c.N_ATTR; i++) {
        if (c.C_1[i].attr == attr)
            g1_copy(c_1, c.C_1[i].c_attr);
    }
}

void get_c_attr_c_2_lu_ok(const uint32_t attr, const struct ciphertext_lu_ok c, g2_t c_2) {
    for (uint32_t i = 0; i < c.N_ATTR; i++) {
        if (c.C_2[i].attr == attr)
            g2_copy(c_2, c.C_2[i].c_attr);
    }
}

void get_c_attr_c_1_ok(const uint32_t attr, const struct ciphertext_ok c, g2_t c_1) {
    for (uint32_t i = 0; i < c.N_ATTR; i++) {
        if (c.C_1[i].attr == attr)
            g2_copy(c_1, c.C_1[i].c_attr);
    }
}

void get_c_attr_c_2_ok(const uint32_t attr, const struct ciphertext_ok c, g2_t c_2) {
    for (uint32_t i = 0; i < c.N_ATTR; i++) {
        if (c.C_2[i].attr == attr)
            g2_copy(c_2, c.C_2[i].c_attr);
    }
}

void get_c_attr_c_1_od(const uint32_t attr, const struct ciphertext_od c, g1_t c_1) {
    for (uint32_t i = 0; i < c.N_ATTR; i++) {
        if (c.C_1[i].attr == attr)
            g1_copy(c_1, c.C_1[i].c_attr);
    }
}

void get_c_attr_c_2_od(const uint32_t attr, const struct ciphertext_od c, g2_t c_2) {
    for (uint32_t i = 0; i < c.N_ATTR; i++) {
        if (c.C_2[i].attr == attr)
            g2_copy(c_2, c.C_2[i].c_attr);
    }
}

void get_c_attr_c_1_cp(const uint32_t attr, const struct ciphertext_cp c, g1_t c_1) {
    for (uint32_t i = 0; i < c.N_ATTR; i++) {
        if (c.C_1[i].attr == attr)
            g1_copy(c_1, c.C_1[i].c_attr);
    }
}

void get_c_attr_c_2_cp(const uint32_t attr, const struct ciphertext_cp c, g2_t c_2) {
    for (uint32_t i = 0; i < c.N_ATTR; i++) {
        if (c.C_2[i].attr == attr)
            g2_copy(c_2, c.C_2[i].c_attr);
    }
}

void set_c_attr_c_1(const uint32_t attr, g1_t c_1, struct ciphertext *c) {
    for (uint32_t i = 0; i < c->N_ATTR; i++) {
        if (c->C_1[i].attr == attr)
            g1_copy(c->C_1[i].c_attr, c_1);
    }
}

void set_c_attr_c_2(const uint32_t attr, g1_t c_2, struct ciphertext *c) {
    for (uint32_t i = 0; i < c->N_ATTR; i++) {
        if (c->C_2[i].attr == attr)
            g1_copy(c->C_2[i].c_attr, c_2);
    }
}

void set_c_attr_c_2_lu_oe(const uint32_t attr, g2_t c_2, struct ciphertext_lu_oe *c) {
    for (uint32_t i = 0; i < c->N_ATTR; i++) {
        if (c->C_2[i].attr == attr)
            g2_copy(c->C_2[i].c_attr, c_2);
    }
}

void set_c_attr_c_1_ok(const uint32_t attr, g2_t c_1, struct ciphertext_ok *c) {
    for (uint32_t i = 0; i < c->N_ATTR; i++) {
        if (c->C_1[i].attr == attr)
            g2_copy(c->C_1[i].c_attr, c_1);
    }
}

void set_c_attr_c_1_lu_ok(const uint32_t attr, g1_t c_1, struct ciphertext_lu_ok *c) {
    for (uint32_t i = 0; i < c->N_ATTR; i++) {
        if (c->C_1[i].attr == attr)
            g1_copy(c->C_1[i].c_attr, c_1);
    }
}

void set_c_attr_c_2_lu_ok(const uint32_t attr, g2_t c_2, struct ciphertext_lu_ok *c) {
    for (uint32_t i = 0; i < c->N_ATTR; i++) {
        if (c->C_2[i].attr == attr)
            g2_copy(c->C_2[i].c_attr, c_2);
    }
}

void set_c_attr_c_2_ok(const uint32_t attr, g2_t c_2, struct ciphertext_ok *c) {
    for (uint32_t i = 0; i < c->N_ATTR; i++) {
        if (c->C_2[i].attr == attr)
            g2_copy(c->C_2[i].c_attr, c_2);
    }
}

void set_c_attr_c_1_od(const uint32_t attr, g1_t c_1, struct ciphertext_od *c) {
    for (uint32_t i = 0; i < c->N_ATTR; i++) {
        if (c->C_1[i].attr == attr)
            g1_copy(c->C_1[i].c_attr, c_1);
    }
}

void set_c_attr_c_2_od(const uint32_t attr, g2_t c_2, struct ciphertext_od *c) {
    for (uint32_t i = 0; i < c->N_ATTR; i++) {
        if (c->C_2[i].attr == attr)
            g2_copy(c->C_2[i].c_attr, c_2);
    }
}

void set_c_attr_c_1_cp(const uint32_t attr, g1_t c_1, struct ciphertext_cp *c) {
    for (uint32_t i = 0; i < c->N_ATTR; i++) {
        if (c->C_1[i].attr == attr)
            g1_copy(c->C_1[i].c_attr, c_1);
    }
}

void set_c_attr_c_1_lu_oe(const uint32_t attr, g1_t c_1, struct ciphertext_lu_oe *c) {
    for (uint32_t i = 0; i < c->N_ATTR; i++) {
        if (c->C_1[i].attr == attr)
            g1_copy(c->C_1[i].c_attr, c_1);
    }
}

void set_c_attr_c_2_cp(const uint32_t attr, g2_t c_2, struct ciphertext_cp *c) {
    for (uint32_t i = 0; i < c->N_ATTR; i++) {
        if (c->C_2[i].attr == attr)
            g2_copy(c->C_2[i].c_attr, c_2);
    }
}

// GPSW
int init_master_key_kp_gpsw(const uint32_t n_attr, struct master_key_kp_gpsw *m) {
    m->N_ATTR = n_attr;
    m->attributes = (struct attribute *)malloc(n_attr * sizeof(struct attribute));

    if (m->attributes == NULL) {
        return EXIT_FAILURE;
    }
    m->t_values = (bn_t *)malloc(n_attr * sizeof(bn_t));
    if (m->t_values == NULL) {
        return EXIT_FAILURE;
    }
    bn_null(m->y);
    bn_new(m->y);
    return EXIT_SUCCESS;
}

int init_public_key_kp_gpsw(const uint32_t n_attr, struct public_key_kp_gpsw *p) {
    p->N_ATTR = n_attr;
    p->T_values = (g2_t *)malloc(n_attr * sizeof(g2_t));
    if (p->T_values == NULL) {
        return EXIT_FAILURE;
    }
    gt_null(&p->Y);
    gt_new(&p->Y);
    return EXIT_SUCCESS;
}

int init_secret_key_kp_gpsw(const uint32_t n_attr, struct secret_key_kp_gpsw *sk) {
    sk->D_values = (g1_t *)malloc(n_attr * sizeof(g1_t));
    if (sk->D_values == NULL) {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int init_ciphertext_kp_gpsw(const uint32_t n_attr, struct ciphertext_kp_gpsw *E) {
    // TODO: Implement to work with complicated access structures. Right now it's just all attributes.
    E->E_values = (g2_t *)malloc(n_attr * sizeof(g2_t));
    if (E->E_values == NULL) {
        return EXIT_FAILURE;
    }
    gt_new(E->E_prime);
    gt_null(E->E_prime);
    return EXIT_SUCCESS;
}

// GPSW Large Universe
int init_master_key_kp_gpsw_lu(const uint32_t n_attr, struct master_key_kp_gpsw_lu *m) {
    bn_new(m->y);
    bn_null(m->y);
    return EXIT_SUCCESS;
}

int init_public_key_kp_gpsw_lu(const uint32_t n_attr, struct public_key_kp_gpsw_lu *p) {
    p->t_values = (g2_t *)malloc((n_attr + 1) * sizeof(g2_t));
    if (p->t_values == NULL) {
        return EXIT_FAILURE;
    }
    g1_new(p->g1);
    g1_null(p->g1);
    g2_new(p->g2);
    g2_null(p->g2);
    return EXIT_SUCCESS;
}

int init_secret_key_kp_gpsw_lu(const uint32_t n_attr, struct secret_key_kp_gpsw_lu *sk) {
    sk->D_values = (g2_t *)malloc(n_attr * sizeof(g2_t));
    if (sk->D_values == NULL) {
        return EXIT_FAILURE;
    } else {
        for (size_t i = 0; i < n_attr; i++) {
            g2_new(sk->D_values[i]);
            g2_null(sk->D_values[i]);
        }
    }
    sk->R_values = (g1_t *)malloc(n_attr * sizeof(g1_t));
    if (sk->R_values == NULL) {
        return EXIT_FAILURE;
    } else {
        for (size_t i = 0; i < n_attr; i++) {
            g1_new(sk->R_values[i]);
            g1_null(sk->R_values[i]);
        }
    }
    return EXIT_SUCCESS;
}

int init_ciphertext_kp_gpsw_lu(const uint32_t n_attr, struct ciphertext_kp_gpsw_lu *E) {
    E->E_values = (g2_t *)malloc(n_attr * sizeof(g2_t));
    if (E->E_values == NULL) {
        return EXIT_FAILURE;
    } else {
        for (size_t i = 0; i < n_attr; i++) {
            g2_new(E->E_values[i]);
            g2_null(E->E_values[i]);
        }
    }
    gt_new(E->E_prime);
    gt_null(E->E_prime);
    g1_new(E->E_prime_prime);
    g1_null(E->E_prime_prime);
    return EXIT_SUCCESS;
}

int init_master_key_k_lin(const uint32_t n_attr, const uint32_t kss, struct master_key_k_lin *m) {
    m->N_ATTR = n_attr;
    m->N_SEC = kss;
    m->atts = (struct k_lin_att * ) malloc ((n_attr + 1) * sizeof(struct k_lin_att));
    m->v_share = (bn_t * ) malloc ((kss+1) * sizeof(bn_t));

    if (m->atts == NULL || m->v_share == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_public_key_k_lin(const uint32_t n_attr, const uint32_t kss, struct public_key_k_lin *p) {
    p->N_ATTR = n_attr;
    p->K_SEC = kss;
    p->mats = (struct k_lin_mat * ) malloc ((2 * (n_attr + 1)) * sizeof(struct k_lin_mat));
    p->a_mat = (g1_t * ) malloc ((kss*(kss+1)) * sizeof(g1_t));
    p->e_mat = (gt_t * ) malloc (kss * sizeof(gt_t));

    if (p->mats == NULL || p->a_mat == NULL || p->e_mat == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_secret_key_K_Lin(const uint32_t n_attr, struct secret_key_K_Lin *s) {
    s->N_ATTR = n_attr;
    s->sk = (struct k_lin_secret_key * ) malloc (((n_attr + 1) * (kss+1)) * sizeof(struct k_lin_secret_key));

    if (s->sk == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_sk_tmp_vj(const uint32_t n_attr, const uint32_t kss, struct sk_tmp_vj *v) {
    v->N_ATTR = n_attr;
    v->K_SEC = kss;
    v->vj = (struct tmp_vj * ) malloc (((n_attr + 1) * (kss+1)) * sizeof(struct tmp_vj));
    v->rj = (struct tmp_rj * ) malloc (((n_attr + 1) * (kss+1)) * sizeof(struct tmp_rj));

    if (v->vj == NULL || v->rj == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_ciphertext_K_Lin(const uint32_t n_attr, const uint32_t kss, struct ciphertext_K_Lin *c) {
    c->N_ATTR = n_attr;
    c->K_SEC = kss;
    c->C_2 = (struct c_attribute_K_Lin * ) malloc ((n_attr + 1) * sizeof(struct c_attribute_K_Lin));
    c->C_1 = (g1_t * ) malloc ((kss+1) * sizeof(g1_t));

    if (c->C_2 == NULL || c->C_1 == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}


//Stuff for klin large universe:
int init_master_key_k_lin_lu(const uint32_t n_attr, const uint32_t kss, struct master_key_k_lin_lu *m) {
    m->N_ATTR = n_attr;
    m->N_SEC = kss;

    m->W_matrix = (bn_t * ) malloc ((((2*kss)+1)*kss) * sizeof(bn_t));
    m->W0_matrix = (bn_t * ) malloc ((((2*kss)+1)*kss) * sizeof(bn_t));
    m->W1_matrix = (bn_t * ) malloc ((((2*kss)+1)*kss) * sizeof(bn_t));
    m->v_secret = (bn_t * ) malloc (((2*kss)+1) * sizeof(bn_t));

    if (m->W_matrix == NULL || m->v_secret == NULL || m->W0_matrix == NULL || m->W1_matrix == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_public_key_k_lin_lu(const uint32_t n_attr, const uint32_t kss, struct public_key_k_lin_lu *p) {
    p->N_ATTR = n_attr;
    p->K_SEC = kss;

    p->A1_mat = (g1_t * ) malloc ((((2*kss)+1)*kss) * sizeof(g1_t));
    p->AW_mat = (g1_t * ) malloc ((kss*kss) * sizeof(g1_t));
    p->AW0_mat = (g1_t * ) malloc ((kss*kss) * sizeof(g1_t));
    p->AW1_mat = (g1_t * ) malloc ((kss*kss) * sizeof(g1_t));
    p->e_mat = (gt_t * ) malloc (kss * sizeof(gt_t));

    if (p->A1_mat == NULL  || p->AW_mat == NULL  || p->AW0_mat == NULL  || p->AW1_mat == NULL  || p->e_mat == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_secret_key_K_Lin_lu(const uint32_t n_attr, struct secret_key_K_Lin_lu *s) {
    s->N_ATTR = n_attr;
    s->sk13 = (struct k_lin_secret_key_lu_13 * ) malloc (((n_attr) * 2*((2*kss)+1) + kss) * sizeof(struct k_lin_secret_key_lu_13));
    s->sk4 = (struct k_lin_secret_key_lu_4 * ) malloc (((n_attr) * 2*((2*kss)+1)) * sizeof(struct k_lin_secret_key_lu_4));

    if (s->sk13 == NULL || s->sk4 == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_sk_tmp_vectors_lu(const uint32_t n_attr, const uint32_t kss, struct sk_tmp_vectors_lu *v) {
    v->N_ATTR = n_attr;
    v->K_SEC = kss;

    v->vj = (struct tmp_vj_lu * ) malloc (((n_attr) * ((2*kss)+1)) * sizeof(struct tmp_vj_lu));
    v->rj = (struct tmp_rj_lu * ) malloc (((n_attr) * ((2*kss)+1)) * sizeof(struct tmp_rj_lu));

    if (v->vj == NULL || v->rj == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_ciphertext_K_Lin_lu(const uint32_t n_attr, const uint32_t kss, struct ciphertext_K_Lin_lu *c) {
    c->N_ATTR = n_attr;
    c->K_SEC = kss;

    c->C_23 = (struct c_attribute_K_Lin_lu_c23 * ) malloc ((n_attr * (((2*kss)+1)+kss)) * sizeof(struct c_attribute_K_Lin_lu_c23));
    c->C_1 = (g1_t * ) malloc (((2*kss)+1) * sizeof(g1_t));

    if (c->C_23 == NULL || c->C_1 == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_tmp_si_lu(const uint32_t n_attr, const uint32_t kss, struct tmp_si_lu *si) {
    si->N_ATTR = n_attr;
    si->K_SEC = kss;

    si->si = (struct si_lu * ) malloc ((n_attr * kss) * sizeof(struct si_lu));

    if (si->si == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}