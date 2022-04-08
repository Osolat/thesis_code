//
// Created by jonas on 2/28/22.
//

#ifndef THESIS_CODE_K_LIN_UTIL_H
#define THESIS_CODE_K_LIN_UTIL_H

#include <stdint.h>

extern "C" {
#include <relic/relic.h>
}

bn_t *vector_trans_mul_matrix(bn_t out[], bn_t v[], bn_t A[], int v_cols, int A_cols, int A_rows, bn_t o);
bn_t *matrix_mul_vector(bn_t out[], bn_t A[], bn_t v[], int a_rows, int a_cols, int v_rows, bn_t o);
bn_t *matrixA_mul_matrixW(bn_t out[], bn_t A[], bn_t Wi[], int a_rows, int a_cols, int w_rows, int w_cols, bn_t o);
bn_t *matrix_mul_scalar(bn_t out[], bn_t A[], int scalar, int a_rows, int a_cols, bn_t o);
bn_t *matrix_add_matrix(bn_t out[], bn_t A[], bn_t Wi[], int a_rows, int a_cols, int w_rows, int w_cols, bn_t o);
bn_t *vector_add_vector(bn_t out[], bn_t v1[], bn_t v2[], int v1_size, int v2_size, bn_t order);
void vector_dot_product(bn_t out, bn_t v1[], bn_t v2[], int v1_size, int v2_size, bn_t o);

g1_t *vector_trans_mul_matrix_g1(g1_t out[], bn_t v[], g1_t A[], int v_cols, int A_cols, int A_rows);
g1_t *matrix_mul_scalar_g1(g1_t out[], g1_t A[], int scalar, int a_rows, int a_cols);
g1_t *vector_add_vector_g1(g1_t out[], g1_t v1[], g1_t v2[], int v1_size, int v2_size);
g1_t *matrix_add_matrix_g1(g1_t out[], g1_t A[], g1_t Wi[], int a_rows, int a_cols, int w_rows, int w_cols);
g1_t *matrixG1_mul_vectorBN(g1_t out[], g1_t A[], bn_t v[], int a_rows, int a_cols, int v_rows);
g1_t *matrixG1_mul_matrixBN(g1_t out[], g1_t A[], bn_t Wi[], int a_rows, int a_cols, int w_rows, int w_cols, g1_t oneValue);

g2_t *vector_trans_mul_matrix_g2(g2_t out[], bn_t v[], g2_t A[], int v_cols, int A_cols, int A_rows);
g2_t *matrix_mul_scalar_g2(g2_t out[], g2_t A[], int scalar, int a_rows, int a_cols);
g2_t *vector_add_vector_g2(g2_t out[], g2_t v1[], g2_t v2[], int v1_size, int v2_size);
g2_t *matrix_add_matrix_g2(g2_t out[], g2_t A[], g2_t Wi[], int a_rows, int a_cols, int w_rows, int w_cols);
g2_t *matrixG2_mul_vectorBN(g2_t out[], g2_t A[], bn_t v[], int a_rows, int a_cols, int v_rows);
g2_t *matrixG2_mul_matrixBN(g2_t out[], g2_t A[], bn_t Wi[], int a_rows, int a_cols, int w_rows, int w_cols, g2_t oneValue);

g1_t *vector_trans_mul_matrix_g1_sim(g1_t out[], bn_t v[], g1_t A[], int v_cols, int A_cols, int A_rows);
g2_t *vector_trans_mul_matrix_g2_sim(g2_t out[], bn_t v[], g2_t A[], int v_cols, int A_cols, int A_rows);

g1_t *vector_trans_mul_matrix_g1_pre(g1_t out[], bn_t v[], g1_t A[][RLC_EP_TABLE_MAX], int v_cols, int A_cols, int A_rows);
g1_t *matrix_mul_scalar_g1_pre(g1_t out[], g1_t A[][RLC_EP_TABLE_MAX], int scalar, int a_rows, int a_cols);
g2_t *vector_trans_mul_matrix_g2_pre(g2_t out[], bn_t v[], g2_t A[][RLC_EP_TABLE_MAX], int v_cols, int A_cols, int A_rows);
g2_t *matrix_mul_scalar_g2_pre(g2_t out[], g2_t A[][RLC_EP_TABLE_MAX], int scalar, int a_rows, int a_cols);

void init_null_new_bn_t_var(bn_t var_name);
void init_null_new_g1_t_var(g1_t var_name);
void init_null_new_g2_t_var(g2_t var_name);
void init_null_new_gt_t_var(gt_t var_name);

void print_sk(struct secret_key_K_Lin *sk, struct sk_tmp_vj *shares, const uint32_t n, const uint32_t ksec);
void print_msk(struct master_key_k_lin *msk, const uint32_t n, const uint32_t ksec);
void print_mpk(struct public_key_k_lin *mpk, const uint32_t n, const uint32_t ksec);
void print_ct(struct ciphertext_K_Lin *CT_A, const uint32_t n, const uint32_t ksec);
void test_matrix_mul_vector(const uint32_t ksec, bn_t o);
void test_vector_trans_mul_matrix(const uint32_t ksec, bn_t o);
void test_matrix_mul_matrix(const uint32_t ksec, bn_t o);
void test_vector_dot_product(const uint32_t ksec, bn_t o);
void test_matrix_mul_scalar(const uint32_t ksec, bn_t o);
void test_matrix_add_matrix(const uint32_t ksec, bn_t o);
void test_vector_add_vector(const uint32_t ksec, bn_t o);
#endif //THESIS_CODE_K_LIN_UTIL_H
