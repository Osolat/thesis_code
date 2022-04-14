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

// GPSW
int init_master_key_kp_gpsw(const uint32_t n_attr, struct master_key_kp_gpsw *m) {
    m->N_ATTR = n_attr;
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

int init_public_key_kp_gpsw_oe(const uint32_t n_attr, struct public_key_kp_gpsw_oe *p) {
    p->N_ATTR = n_attr;
    p->T_values = (g1_t *)malloc(n_attr * sizeof(g1_t));
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

int init_secret_key_kp_gpsw_oe(const uint32_t n_attr, struct secret_key_kp_gpsw_oe *sk) {
    sk->D_values = (g2_t *)malloc(n_attr * sizeof(g2_t));
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

int init_ciphertext_kp_gpsw_oe(const uint32_t n_attr, struct ciphertext_kp_gpsw_oe *E) {
    // TODO: Implement to work with complicated access structures. Right now it's just all attributes.
    E->E_values = (g1_t *)malloc(n_attr * sizeof(g1_t));
    if (E->E_values == NULL) {
        return EXIT_FAILURE;
    }
    gt_new(E->E_prime);
    gt_null(E->E_prime);
    return EXIT_SUCCESS;
}

// GPSW Large Universe
int init_master_key_kp_gpsw_lu_ok(const uint32_t n_attr, struct master_key_kp_gpsw_lu_ok *m) {
    bn_new(m->y);
    bn_null(m->y);
    return EXIT_SUCCESS;
}

int init_public_key_kp_gpsw_lu_ok(const uint32_t n_attr, struct public_key_kp_gpsw_lu_ok *p) {
    p->t_values = (g2_t *)malloc((n_attr + 1) * sizeof(g2_t));
    if (p->t_values == NULL) {
        return EXIT_FAILURE;
    } else {
        for (size_t i = 0; i < n_attr + 1; i++) {
            g2_new(p->t_values[i]);
            g2_null(p->t_values[i]);
        }
    }
    g1_new(p->g1);
    g1_null(p->g1);
    g2_new(p->g2);
    g2_null(p->g2);
    return EXIT_SUCCESS;
}

int init_secret_key_kp_gpsw_lu_ok(const uint32_t n_attr, struct secret_key_kp_gpsw_lu_ok *sk) {
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

int init_ciphertext_kp_gpsw_lu_ok(const uint32_t n_attr, struct ciphertext_kp_gpsw_lu_ok *E) {
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
    m->atts = (struct k_lin_att *)malloc((n_attr + 1) * sizeof(struct k_lin_att));
    m->v_share = (bn_t *)malloc((kss + 1) * sizeof(bn_t));

    if (m->atts == NULL || m->v_share == NULL) {
        return EXIT_FAILURE;
    } else {
        for (size_t i = 0; i < n_attr; i++) {
            for (int j = 0; j < ((kss + 1) * kss); ++j) {
                bn_null(m->atts[i].w[j]);
                bn_new(m->atts[i].w[j]);
            }
        }
        return EXIT_SUCCESS;
    }
}

int init_public_key_k_lin(const uint32_t n_attr, const uint32_t kss, struct public_key_k_lin *p) {
    p->N_ATTR = n_attr;
    p->K_SEC = kss;
    p->mats = (struct k_lin_mat *)malloc((2 * (n_attr + 1)) * sizeof(struct k_lin_mat));
    p->a_mat = (g1_t *)malloc((kss * (kss + 1)) * sizeof(g1_t));
    p->e_mat = (gt_t *)malloc(kss * sizeof(gt_t));

    if (p->mats == NULL || p->a_mat == NULL || p->e_mat == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_secret_key_K_Lin(const uint32_t n_attr, struct secret_key_K_Lin *s) {
    s->N_ATTR = n_attr;
    s->sk = (struct k_lin_secret_key *)malloc(((n_attr + 1) * (kss + 1)) * sizeof(struct k_lin_secret_key));

    if (s->sk == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_sk_tmp_vj(const uint32_t n_attr, const uint32_t kss, struct sk_tmp_vj *v) {
    v->N_ATTR = n_attr;
    v->K_SEC = kss;
    v->vj = (struct tmp_vj *)malloc(((n_attr + 1) * (kss + 1)) * sizeof(struct tmp_vj));
    v->rj = (struct tmp_rj *)malloc(((n_attr + 1) * (kss + 1)) * sizeof(struct tmp_rj));

    if (v->vj == NULL || v->rj == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_ciphertext_K_Lin(const uint32_t n_attr, const uint32_t kss, struct ciphertext_K_Lin *c) {
    c->N_ATTR = n_attr;
    c->K_SEC = kss;
    c->C_2 = (struct c_attribute_K_Lin *)malloc((n_attr + 1) * sizeof(struct c_attribute_K_Lin));
    c->C_1 = (g1_t *)malloc((kss + 1) * sizeof(g1_t));

    if (c->C_2 == NULL || c->C_1 == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

// Stuff for klin large universe:
int init_master_key_k_lin_lu(const uint32_t n_attr, const uint32_t kss, struct master_key_k_lin_lu *m) {
    m->N_ATTR = n_attr;
    m->N_SEC = kss;

    m->W_matrix = (bn_t *)malloc((((2 * kss) + 1) * kss) * sizeof(bn_t));
    m->W0_matrix = (bn_t *)malloc((((2 * kss) + 1) * kss) * sizeof(bn_t));
    m->W1_matrix = (bn_t *)malloc((((2 * kss) + 1) * kss) * sizeof(bn_t));
    m->v_secret = (bn_t *)malloc(((2 * kss) + 1) * sizeof(bn_t));

    if (m->W_matrix == NULL || m->v_secret == NULL || m->W0_matrix == NULL || m->W1_matrix == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_public_key_k_lin_lu(const uint32_t n_attr, const uint32_t kss, struct public_key_k_lin_lu *p) {
    p->N_ATTR = n_attr;
    p->K_SEC = kss;

    p->A1_mat = (g1_t *)malloc((((2 * kss) + 1) * kss) * sizeof(g1_t));
    p->AW_mat = (g1_t *)malloc((kss * kss) * sizeof(g1_t));
    p->AW0_mat = (g1_t *)malloc((kss * kss) * sizeof(g1_t));
    p->AW1_mat = (g1_t *)malloc((kss * kss) * sizeof(g1_t));
    p->e_mat = (gt_t *)malloc(kss * sizeof(gt_t));

    if (p->A1_mat == NULL || p->AW_mat == NULL || p->AW0_mat == NULL || p->AW1_mat == NULL || p->e_mat == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_secret_key_K_Lin_lu(const uint32_t n_attr, struct secret_key_K_Lin_lu *s) {
    s->N_ATTR = n_attr;
    s->sk13 = (struct k_lin_secret_key_lu_13 *)malloc(((n_attr)*2 * ((2 * kss) + 1) + kss) * sizeof(struct k_lin_secret_key_lu_13));
    s->sk4 = (struct k_lin_secret_key_lu_4 *)malloc(((n_attr)*2 * ((2 * kss) + 1)) * sizeof(struct k_lin_secret_key_lu_4));

    if (s->sk13 == NULL || s->sk4 == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_sk_tmp_vectors_lu(const uint32_t n_attr, const uint32_t kss, struct sk_tmp_vectors_lu *v) {
    v->N_ATTR = n_attr;
    v->K_SEC = kss;

    v->vj = (struct tmp_vj_lu *)malloc(((n_attr) * ((2 * kss) + 1)) * sizeof(struct tmp_vj_lu));
    v->rj = (struct tmp_rj_lu *)malloc(((n_attr) * ((2 * kss) + 1)) * sizeof(struct tmp_rj_lu));

    if (v->vj == NULL || v->rj == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_ciphertext_K_Lin_lu(const uint32_t n_attr, const uint32_t kss, struct ciphertext_K_Lin_lu *c) {
    c->N_ATTR = n_attr;
    c->K_SEC = kss;

    c->C_23 = (struct c_attribute_K_Lin_lu_c23 *)malloc((n_attr * (((2 * kss) + 1) + kss)) * sizeof(struct c_attribute_K_Lin_lu_c23));
    c->C_1 = (g1_t *)malloc(((2 * kss) + 1) * sizeof(g1_t));

    if (c->C_23 == NULL || c->C_1 == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_tmp_si_lu(const uint32_t n_attr, const uint32_t kss, struct tmp_si_lu *si) {
    si->N_ATTR = n_attr;
    si->K_SEC = kss;

    si->si = (struct si_lu *)malloc((n_attr * kss) * sizeof(struct si_lu));

    if (si->si == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

// Optimised KeyGen
int init_master_key_k_lin_ok(const uint32_t n_attr, const uint32_t kss, struct master_key_k_lin_ok *m) {
    m->N_ATTR = n_attr;
    m->N_SEC = kss;
    m->atts = (struct k_lin_att_ok *)malloc((n_attr + 1) * sizeof(struct k_lin_att_ok));
    m->v_share = (bn_t *)malloc((kss + 1) * sizeof(bn_t));

    if (m->atts == NULL || m->v_share == NULL) {
        return EXIT_FAILURE;
    } else {
        for (size_t i = 0; i < n_attr; i++) {
            for (int j = 0; j < ((kss + 1) * kss); ++j) {
                bn_null(m->atts[i].w[j]);
                bn_new(m->atts[i].w[j]);
            }
        }
        return EXIT_SUCCESS;
    }
}

int init_public_key_k_lin_ok(const uint32_t n_attr, const uint32_t kss, struct public_key_k_lin_ok *p) {
    p->N_ATTR = n_attr;
    p->K_SEC = kss;
    p->mats = (struct k_lin_mat_ok *)malloc((2 * (n_attr + 1)) * sizeof(struct k_lin_mat_ok));
    p->a_mat = (g2_t *)malloc((kss * (kss + 1)) * sizeof(g2_t));
    p->e_mat = (gt_t *)malloc(kss * sizeof(gt_t));

    if (p->mats == NULL || p->a_mat == NULL || p->e_mat == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_secret_key_K_Lin_ok(const uint32_t n_attr, struct secret_key_K_Lin_ok *s) {
    s->N_ATTR = n_attr;
    s->sk = (struct k_lin_secret_key_ok *)malloc(((n_attr + 1) * (kss + 1)) * sizeof(struct k_lin_secret_key_ok));

    if (s->sk == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_sk_tmp_vj_ok(const uint32_t n_attr, const uint32_t kss, struct sk_tmp_vj_ok *v) {
    v->N_ATTR = n_attr;
    v->K_SEC = kss;
    v->vj = (struct tmp_vj_ok *)malloc(((n_attr + 1) * (kss + 1)) * sizeof(struct tmp_vj_ok));
    v->rj = (struct tmp_rj_ok *)malloc(((n_attr + 1) * (kss + 1)) * sizeof(struct tmp_rj_ok));

    if (v->vj == NULL || v->rj == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_ciphertext_K_Lin_ok(const uint32_t n_attr, const uint32_t kss, struct ciphertext_K_Lin_ok *c) {
    c->N_ATTR = n_attr;
    c->K_SEC = kss;
    c->C_2 = (struct c_attribute_K_Lin_ok *)malloc((n_attr + 1) * sizeof(struct c_attribute_K_Lin_ok));
    c->C_1 = (g2_t *)malloc((kss + 1) * sizeof(g2_t));

    if (c->C_2 == NULL || c->C_1 == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

// Stuff for optimised decryption
// TODO move stuff around in groups, move ct1 to h and sk1 to g. Other changes may be needed.
int init_master_key_k_lin_od(const uint32_t n_attr, const uint32_t kss, struct master_key_k_lin_od *m) {
    m->N_ATTR = n_attr;
    m->N_SEC = kss;
    m->atts = (struct k_lin_att_od *)malloc((n_attr + 1) * sizeof(struct k_lin_att_od));
    m->v_share = (bn_t *)malloc((kss + 1) * sizeof(bn_t));

    if (m->atts == NULL || m->v_share == NULL) {
        return EXIT_FAILURE;
    } else {
        for (size_t i = 0; i < n_attr; i++) {
            for (int j = 0; j < ((kss + 1) * kss); ++j) {
                bn_null(m->atts[i].w[j]);
                bn_new(m->atts[i].w[j]);
            }
        }
        return EXIT_SUCCESS;
    }
}

int init_public_key_k_lin_od(const uint32_t n_attr, const uint32_t kss, struct public_key_k_lin_od *p) {
    p->N_ATTR = n_attr;
    p->K_SEC = kss;
    p->mats = (struct k_lin_mat_od *)malloc((2 * (n_attr + 1)) * sizeof(struct k_lin_mat_od));
    p->a_mat = (g2_t *)malloc((kss * (kss + 1)) * sizeof(g2_t));
    p->e_mat = (gt_t *)malloc(kss * sizeof(gt_t));

    if (p->mats == NULL || p->a_mat == NULL || p->e_mat == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_secret_key_K_Lin_od(const uint32_t n_attr, struct secret_key_K_Lin_od *s) {
    s->N_ATTR = n_attr;
    s->sk = (struct k_lin_secret_key_od *)malloc(((n_attr + 1) * (kss + 1)) * sizeof(struct k_lin_secret_key_od));

    if (s->sk == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_sk_tmp_vj_od(const uint32_t n_attr, const uint32_t kss, struct sk_tmp_vj_od *v) {
    v->N_ATTR = n_attr;
    v->K_SEC = kss;
    v->vj = (struct tmp_vj_od *)malloc(((n_attr + 1) * (kss + 1)) * sizeof(struct tmp_vj_od));
    v->rj = (struct tmp_rj_od *)malloc(((n_attr + 1) * (kss + 1)) * sizeof(struct tmp_rj_od));

    if (v->vj == NULL || v->rj == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_ciphertext_K_Lin_od(const uint32_t n_attr, const uint32_t kss, struct ciphertext_K_Lin_od *c) {
    c->N_ATTR = n_attr;
    c->K_SEC = kss;
    c->C_2 = (struct c_attribute_K_Lin_od *)malloc((n_attr + 1) * sizeof(struct c_attribute_K_Lin_od));
    c->C_1 = (g2_t *)malloc((kss + 1) * sizeof(g2_t));

    if (c->C_2 == NULL || c->C_1 == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}


// Stuff for klin unbounded version OD:
int init_master_key_k_lin_ub_od(const uint32_t n_attr, const uint32_t kss, struct master_key_k_lin_ub_od *m) {
    m->N_ATTR = n_attr;
    m->N_SEC = kss;

    m->W_matrix = (bn_t *)malloc((((2 * kss) + 1) * kss) * sizeof(bn_t));
    m->W0_matrix = (bn_t *)malloc((((2 * kss) + 1) * kss) * sizeof(bn_t));
    m->W1_matrix = (bn_t *)malloc((((2 * kss) + 1) * kss) * sizeof(bn_t));
    m->v_secret = (bn_t *)malloc(((2 * kss) + 1) * sizeof(bn_t));

    if (m->W_matrix == NULL || m->v_secret == NULL || m->W0_matrix == NULL || m->W1_matrix == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_public_key_k_lin_ub_od(const uint32_t n_attr, const uint32_t kss, struct public_key_k_lin_ub_od *p) {
    p->N_ATTR = n_attr;
    p->K_SEC = kss;

    p->A1_mat_g1 = (g1_t *)malloc((((2 * kss) + 1) * kss) * sizeof(g1_t));
    p->A1_mat_g2 = (g2_t *)malloc((((2 * kss) + 1) * kss) * sizeof(g2_t));
    p->AW_mat = (g1_t *)malloc((kss * kss) * sizeof(g1_t));
    p->AW0_mat = (g1_t *)malloc((kss * kss) * sizeof(g1_t));
    p->AW1_mat = (g1_t *)malloc((kss * kss) * sizeof(g1_t));
    p->e_mat = (gt_t *)malloc(kss * sizeof(gt_t));

    if (p->A1_mat_g1 == NULL || p->A1_mat_g1 == NULL || p->AW_mat == NULL || p->AW0_mat == NULL || p->AW1_mat == NULL || p->e_mat == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_secret_key_K_Lin_ub_od(const uint32_t n_attr, struct secret_key_K_Lin_ub_od *s) {
    s->N_ATTR = n_attr;
    s->sk13 = (struct k_lin_secret_key_ub_13_od *)malloc(((n_attr)*2 * ((2 * kss) + 1) + kss) * sizeof(struct k_lin_secret_key_ub_13_od));
    s->sk4 = (struct k_lin_secret_key_ub_4_od *)malloc(((n_attr)*2 * ((2 * kss) + 1)) * sizeof(struct k_lin_secret_key_ub_4_od));

    if (s->sk13 == NULL || s->sk4 == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_sk_tmp_vectors_ub_od(const uint32_t n_attr, const uint32_t kss, struct sk_tmp_vectors_ub_od *v) {
    v->N_ATTR = n_attr;
    v->K_SEC = kss;

    v->vj = (struct tmp_vj_ub_od *)malloc(((n_attr) * ((2 * kss) + 1)) * sizeof(struct tmp_vj_ub_od));
    v->rj = (struct tmp_rj_ub_od *)malloc(((n_attr) * ((2 * kss) + 1)) * sizeof(struct tmp_rj_ub_od));

    if (v->vj == NULL || v->rj == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_ciphertext_K_Lin_ub_od(const uint32_t n_attr, const uint32_t kss, struct ciphertext_K_Lin_ub_od *c) {
    c->N_ATTR = n_attr;
    c->K_SEC = kss;

    c->C_23 = (struct c_attribute_K_Lin_ub_c23_od *)malloc((n_attr * (((2 * kss) + 1) + kss)) * sizeof(struct c_attribute_K_Lin_ub_c23_od));
    c->C_1 = (g2_t *)malloc(((2 * kss) + 1) * sizeof(g2_t));

    if (c->C_23 == NULL || c->C_1 == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_tmp_si_ub_od(const uint32_t n_attr, const uint32_t kss, struct tmp_si_ub_od *si) {
    si->N_ATTR = n_attr;
    si->K_SEC = kss;

    si->si = (struct si_ub_od *)malloc((n_attr * kss) * sizeof(struct si_ub_od));

    if (si->si == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

// Stuff for klin unbounded version OK:
int init_master_key_k_lin_ub_ok(const uint32_t n_attr, const uint32_t kss, struct master_key_k_lin_ub_ok *m) {
    m->N_ATTR = n_attr;
    m->N_SEC = kss;

    m->W_matrix = (bn_t *)malloc((((2 * kss) + 1) * kss) * sizeof(bn_t));
    m->W0_matrix = (bn_t *)malloc((((2 * kss) + 1) * kss) * sizeof(bn_t));
    m->W1_matrix = (bn_t *)malloc((((2 * kss) + 1) * kss) * sizeof(bn_t));
    m->v_secret = (bn_t *)malloc(((2 * kss) + 1) * sizeof(bn_t));

    if (m->W_matrix == NULL || m->v_secret == NULL || m->W0_matrix == NULL || m->W1_matrix == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_public_key_k_lin_ub_ok(const uint32_t n_attr, const uint32_t kss, struct public_key_k_lin_ub_ok *p) {
    p->N_ATTR = n_attr;
    p->K_SEC = kss;

    p->A1_mat = (g2_t *)malloc((((2 * kss) + 1) * kss) * sizeof(g2_t));
    p->AW_mat = (g2_t *)malloc((kss * kss) * sizeof(g2_t));
    p->AW0_mat = (g2_t *)malloc((kss * kss) * sizeof(g2_t));
    p->AW1_mat = (g2_t *)malloc((kss * kss) * sizeof(g2_t));
    p->e_mat = (gt_t *)malloc(kss * sizeof(gt_t));

    if (p->A1_mat == NULL || p->AW_mat == NULL || p->AW0_mat == NULL || p->AW1_mat == NULL || p->e_mat == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_secret_key_K_Lin_ub_ok(const uint32_t n_attr, struct secret_key_K_Lin_ub_ok *s) {
    s->N_ATTR = n_attr;
    s->sk13 = (struct k_lin_secret_key_ub_13_ok *)malloc(((n_attr)*2 * ((2 * kss) + 1) + kss) * sizeof(struct k_lin_secret_key_ub_13_ok));
    s->sk4 = (struct k_lin_secret_key_ub_4_ok *)malloc(((n_attr)*2 * ((2 * kss) + 1)) * sizeof(struct k_lin_secret_key_ub_4_ok));

    if (s->sk13 == NULL || s->sk4 == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_sk_tmp_vectors_ub_ok(const uint32_t n_attr, const uint32_t kss, struct sk_tmp_vectors_ub_ok *v) {
    v->N_ATTR = n_attr;
    v->K_SEC = kss;

    v->vj = (struct tmp_vj_ub_ok *)malloc(((n_attr) * ((2 * kss) + 1)) * sizeof(struct tmp_vj_ub_ok));
    v->rj = (struct tmp_rj_ub_ok *)malloc(((n_attr) * ((2 * kss) + 1)) * sizeof(struct tmp_rj_ub_ok));

    if (v->vj == NULL || v->rj == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_ciphertext_K_Lin_ub_ok(const uint32_t n_attr, const uint32_t kss, struct ciphertext_K_Lin_ub_ok *c) {
    c->N_ATTR = n_attr;
    c->K_SEC = kss;

    c->C_23 = (struct c_attribute_K_Lin_ub_c23_ok *)malloc((n_attr * (((2 * kss) + 1) + kss)) * sizeof(struct c_attribute_K_Lin_ub_c23_ok));
    c->C_1 = (g2_t *)malloc(((2 * kss) + 1) * sizeof(g2_t));

    if (c->C_23 == NULL || c->C_1 == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_tmp_si_ub_ok(const uint32_t n_attr, const uint32_t kss, struct tmp_si_ub_ok *si) {
    si->N_ATTR = n_attr;
    si->K_SEC = kss;

    si->si = (struct si_ub_ok *)malloc((n_attr * kss) * sizeof(struct si_ub_ok));

    if (si->si == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

//ALP


int init_public_params_alp_oe(const uint32_t bound, struct alp_pp_naive_oe *pp){
    pp->bound = bound;
    pp->U1 = (g1_t *)malloc((bound+1)*sizeof(g1_t));
    pp->U2 = (g2_t *)malloc((bound+1)*sizeof(g2_t));
    pp->H1 = (g1_t *)malloc(bound*sizeof(g1_t));
    pp->H2 = (g2_t *)malloc(bound*sizeof(g2_t));
    g1_null(pp->g1); g1_new(pp->g1); 
    g1_null(pp->g2); g1_null(pp->g2); 
    gt_null(pp->gt); gt_null(pp->gt); 
    if (pp->U1 == NULL || pp->H1 == NULL || pp->H2 == NULL || pp->U2 == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_secret_key_attr_alp_oe(const uint32_t bound, struct alp_sk_attr_oe *sk) { 
    g2_null(sk -> D1); g2_new(sk -> D1); 
    g2_null(sk -> D2); g2_new(sk -> D2);
    uint32_t k_length = bound-1;
    sk -> K = (g2_t *)malloc( k_length*sizeof(g2_t) ); 
    if (sk -> D1 == NULL || sk -> D2 == NULL || sk -> K == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_secret_key_alp_oe(const uint32_t bound, struct alp_sk_oe *sk) {
    sk -> D = (struct alp_sk_attr_oe *)malloc(bound * sizeof(struct alp_sk_attr_oe));
    if (sk -> D == NULL){
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}
