//
// Created by jonas on 2/28/22.
//

#include "k_lin_util.h"
#include <stdlib.h>

#include "../structures.h"
#include "../zp_arith.h"
#include "../g1_arith.h"
#include "../g2_arith.h"
#include "../gt_arith.h"

void init_null_new_bn_t_var(bn_t var_name){
    bn_null(var_name);
    bn_new(var_name);
}

void init_null_new_g1_t_var(g1_t var_name){
    g1_null(var_name);
    g1_new(var_name);
}

void init_null_new_g2_t_var(g2_t var_name){
    g2_null(var_name);
    g2_new(var_name);
}

void init_null_new_gt_t_var(gt_t var_name){
    gt_null(var_name);
    gt_new(var_name);
}

bn_t *vector_trans_mul_matrix(bn_t out[], bn_t v[], bn_t A[], int v_cols, int A_cols, int A_rows, bn_t ord) {
    if (v_cols != A_rows) {
        printf("Error in matrix - trans vector dimensions! \n");
        exit(-1);
    }

    bn_t tmp; bn_t tmp_add;
    init_null_new_bn_t_var(tmp);
    init_null_new_bn_t_var(tmp_add);

    int ctr = 0;
    int z = 0;
    int h = 0;
    for (int i = 0; i < A_cols; ++i) {
        z = h;
        for (int j = 0; j < v_cols; ++j, z += A_cols) {
            bn_t_mul(tmp, v[j], A[z], ord);
            bn_t_add(tmp_add, tmp_add, tmp, ord);
        }
        h++;
        bn_copy(out[ctr], tmp_add);
        bn_t_setzero(tmp_add);
        ctr ++;
    }
    return out;
}

bn_t *matrix_mul_vector(bn_t out[], bn_t A[], bn_t v[], int a_rows, int a_cols, int v_rows, bn_t o) {
    if (a_cols != v_rows) {
        printf("Error in matrix - vector dimensions! \n");
        exit(-1);
    }

    bn_t tmp; bn_t tmp_add;
    init_null_new_bn_t_var(tmp);
    init_null_new_bn_t_var(tmp_add);
    int z = 0;
    int cter = 0;
    int h = 0;
    for (int x = 0; x < (a_rows); x++, h += a_cols) {
        z = h;
        for (int y = 0; y < v_rows; ++y) {
            bn_t_mul(tmp, A[z], v[y], o);
            bn_t_add(tmp_add, tmp_add, tmp, o);
            z++;
        }

        bn_copy(out[cter], tmp_add);
        bn_t_setzero(tmp_add);
        cter ++;
    }
    return out;
}

bn_t *matrixA_mul_matrixW(bn_t out[], bn_t A[], bn_t Wi[], int a_rows, int a_cols, int w_rows, int w_cols, bn_t o) {
    if (a_cols != w_rows) {
        printf("Error in matrix - matrix dimensions! \n");
        exit(-1);
    }

    bn_t tmp; bn_t tmp_add;
    init_null_new_bn_t_var(tmp);
    init_null_new_bn_t_var(tmp_add);

    int ctt = 0;
    int z = 0;
    int h = 0;
    int k = 0;
    int g = 0;
    for (int x = 0; x < a_rows; x++, g += a_cols) {
        //2x
        h = 0;
        for (int i = 0; i < w_cols; ++i) {
            //2x
            z = h;
            k = g;
            for (int j = 0; j < a_cols; ++j, z += w_cols) {
                //x3
                bn_t_mul(tmp, A[k], Wi[z], o);
                bn_t_add(tmp_add, tmp_add, tmp, o);
                k++;
            }
            h++;
            bn_copy(out[ctt], tmp_add);
            bn_t_setzero(tmp_add);
            ctt ++;
        }
        bn_set_dig(tmp, 1);
        bn_t_setzero(tmp_add);
    }
    return out;
}

bn_t *matrix_mul_scalar(bn_t out[], bn_t A[], int scalar, int a_rows, int a_cols, bn_t o) {
    bn_t tmp_mul;
    init_null_new_bn_t_var(tmp_mul);
    bn_t scalarBn_t;
    init_null_new_bn_t_var(scalarBn_t);
    bn_set_dig(scalarBn_t, scalar);

    for (int i = 0; i < (a_rows * a_cols); ++i) {
        bn_t_mul(tmp_mul, A[i], scalarBn_t, o);
        bn_copy(out[i], tmp_mul);
        bn_t_setzero(tmp_mul);
    }
    return out;
}

bn_t *matrix_add_matrix(bn_t out[], bn_t A[], bn_t Wi[], int a_rows, int a_cols, int w_rows, int w_cols, bn_t o) {
    if (a_rows != w_rows || a_cols != w_cols) {
        printf("Error in matrix add matrix dimensions! \n");
        exit(-1);
    }
    bn_t tmp_add;
    init_null_new_bn_t_var(tmp_add);

    for (int i = 0; i < (a_rows * a_cols); ++i) {
        bn_t_add(tmp_add, A[i], Wi[i], o);
        bn_copy(out[i], tmp_add);
        bn_t_setzero(tmp_add);
    }

    return out;
}

bn_t *vector_add_vector(bn_t out[], bn_t v1[], bn_t v2[], int v1_size, int v2_size, bn_t order) {
    if (v1_size != v2_size) {
        printf("Error in vector add vector dimensions! \n");
        exit(-1);
    }
    bn_t tmp_add;
    init_null_new_bn_t_var(tmp_add);
    for (int i = 0; i < v1_size; ++i) {
        bn_t_add(tmp_add, v1[i], v2[i], order);
        bn_copy(out[i], tmp_add);
        bn_t_setzero(tmp_add);
        //bn_null(tmp_add); bn_new(tmp_add);
    }
    return out;
}

void vector_dot_product(bn_t out, bn_t v1[], bn_t v2[], int v1_size, int v2_size, bn_t o) {
    if (v1_size != v2_size) {
        printf("Error in vector - dot - vector dimensions! \n");
        exit(-1);
    }
    bn_t tmp;
    bn_t tmp_add;
    init_null_new_bn_t_var(tmp);
    init_null_new_bn_t_var(tmp_add);
    for (int i = 0; i < v1_size; ++i) {
        bn_t_mul(tmp, v1[i], v2[i], o);
        bn_t_add(tmp_add, tmp_add, tmp, o);
    }
    bn_copy(out,tmp_add);
}




//g1 versions
g1_t *vector_trans_mul_matrix_g1(g1_t out[], bn_t v[], g1_t A[], int v_cols, int A_cols, int A_rows) {
    if (v_cols != A_rows) {
        printf("Error in matrix - trans vector dimensions! \n");
        exit(-1);
    }

    g1_t tmp; g1_t tmp_add;
    init_null_new_g1_t_var(tmp);
    init_null_new_g1_t_var(tmp_add);
    g1_mul_dig(tmp_add, tmp_add, 0);

    int ctr = 0;
    int z = 0;
    int h = 0;
    for (int i = 0; i < A_cols; ++i) {
        z = h;
        for (int j = 0; j < v_cols; ++j, z += A_cols) {
            g1_mul(tmp, A[z], v[j]);
            g1_add(tmp_add, tmp_add, tmp);
        }
        h++;
        g1_copy(out[ctr], tmp_add);
        g1_mul_dig(tmp_add, tmp_add, 0);
        ctr ++;
    }
    return out;
}

g1_t *matrix_mul_scalar_g1(g1_t out[], g1_t A[], int scalar, int a_rows, int a_cols) {
    g1_t tmp_mul;
    init_null_new_g1_t_var(tmp_mul);
    g1_mul_dig(tmp_mul, tmp_mul, 0);

    bn_t scalarBn_t;
    init_null_new_bn_t_var(scalarBn_t);
    bn_set_dig(scalarBn_t, scalar);

    for (int i = 0; i < (a_rows * a_cols); ++i) {
        g1_mul(tmp_mul, A[i], scalarBn_t);
        g1_copy(out[i], tmp_mul);
        g1_mul_dig(tmp_mul, tmp_mul, 0);
    }
    return out;
}

g1_t *vector_add_vector_g1(g1_t out[], g1_t v1[], g1_t v2[], int v1_size, int v2_size) {
    if (v1_size != v2_size) {
        printf("Error in vector add vector dimensions! \n");
        exit(-1);
    }
    g1_t tmp_add;
    init_null_new_g1_t_var(tmp_add);
    g1_mul_dig(tmp_add, tmp_add, 0);

    for (int i = 0; i < v1_size; ++i) {
        g1_add(tmp_add, v1[i], v2[i]);
        g1_copy(out[i], tmp_add);
        g1_mul_dig(tmp_add, tmp_add, 0);
    }
    return out;
}

g1_t *matrix_add_matrix_g1(g1_t out[], g1_t A[], g1_t Wi[], int a_rows, int a_cols, int w_rows, int w_cols) {
    if (a_rows != w_rows || a_cols != w_cols) {
        printf("Error in matrix add matrix dimensions! \n");
        exit(-1);
    }
    g1_t tmp_add;
    init_null_new_g1_t_var(tmp_add);
    g1_mul_dig(tmp_add, tmp_add, 0);

    for (int i = 0; i < (a_rows * a_cols); ++i) {
        g1_add(tmp_add, A[i], Wi[i]);
        g1_copy(out[i], tmp_add);
        g1_mul_dig(tmp_add, tmp_add, 0);
    }

    return out;
}

g1_t *matrixG1_mul_vectorBN(g1_t out[], g1_t A[], bn_t v[], int a_rows, int a_cols, int v_rows) {
    if (a_cols != v_rows) {
        printf("Error in matrix - vector dimensions! \n");
        exit(-1);
    }

    g1_t tmp; g1_t tmp_add;
    init_null_new_g1_t_var(tmp);
    init_null_new_g1_t_var(tmp_add);
    g1_mul_dig(tmp_add, tmp_add, 0);

    int z = 0;
    int cter = 0;
    int h = 0;
    for (int x = 0; x < (a_rows); x++, h += a_cols) {
        z = h;
        for (int y = 0; y < v_rows; ++y) {
            g1_mul(tmp, A[z], v[y]);
            g1_add(tmp_add, tmp_add, tmp);
            z++;
        }

        g1_copy(out[cter], tmp_add);
        g1_mul_dig(tmp_add, tmp_add, 0);
        cter ++;
    }
    return out;
}

g1_t *matrixG1_mul_matrixBN(g1_t out[], g1_t A[], bn_t Wi[], int a_rows, int a_cols, int w_rows, int w_cols, g1_t oneValue) {
    if (a_cols != w_rows) {
        printf("Error in matrix - matrix dimensions! \n");
        exit(-1);
    }

    g1_t tmp; g1_t tmp_add;
    init_null_new_g1_t_var(tmp);
    init_null_new_g1_t_var(tmp_add);
    g1_mul_dig(tmp_add, tmp_add, 0);

    int ctt = 0;
    int z = 0;
    int h = 0;
    int k = 0;
    int g = 0;
    for (int x = 0; x < a_rows; x++, g += a_cols) {
        //2x
        h = 0;
        for (int i = 0; i < w_cols; ++i) {
            //2x
            z = h;
            k = g;
            for (int j = 0; j < a_cols; ++j, z += w_cols) {
                //x3
                g1_mul(tmp, A[k], Wi[z]);
                g1_add(tmp_add, tmp_add, tmp);
                k++;
            }
            h++;
            g1_copy(out[ctt], tmp_add);
            g1_mul_dig(tmp_add, tmp_add, 0);
            ctt ++;
        }
        g1_mul_dig(tmp, tmp, 0);
        g1_add(tmp, tmp, oneValue);
        g1_mul_dig(tmp_add, tmp_add, 0);
    }
    return out;
}


void print_sk(struct secret_key_K_Lin *sk, struct sk_tmp_vj *shares, const uint32_t n, const uint32_t ksec) {
    for (int j = 0; j < (n + 1); ++j) {
        printf("sk_1_%d: \n", j);
        for (int i = 0; i < (ksec + 1); ++i) {
            printf("sk_1_%d[%d]: \n", j, i);
            g2_print(sk->sk[j].sk_one[i]);
        }
    }

    printf("\n");

    for (int x = 0; x < (n + 1); ++x) {
        printf("sk_2_%d: \n", x);
        for (int y = 0; y < ksec; ++y) {
            printf("sk_2_%d[%d]: \n", x, y);
            g2_print(sk->sk[x].sk_two[y]);
        }
    }

    printf("\n");

    for (int a = 0; a < n; ++a) {
        printf("vector vj_%d: \n", a);
        for (int b = 0; b < (ksec + 1); ++b) {
            printf("vector vj_%d[%d]: ", a, b);
            bn_print(shares->vj[a].vec_j[b]);
        }
    }

    printf("\n");

    for (int c = 0; c < n; ++c) {
        printf("vector rj_%d: \n", c);
        for (int d = 0; d < (ksec); ++d) {
            printf("vector rj_%d[%d]: ", c, d);
            bn_print(shares->rj[c].vec_rj[d]);
        }
    }


}

void print_msk(struct master_key_k_lin *msk, const uint32_t n, const uint32_t ksec) {
    for (int i = 0; i < (ksec + 1); ++i) {
        printf("Vector v[%d]: ", i);
        bn_print(msk->v_share[i]);
    }
    printf("\n");
    for (int j = 0; j < (n+1); ++j) {
        printf("Matrix W_%d: \n", j);
        for (int x = 0; x < ((ksec + 1) * ksec); ++x) {
            printf("Matrix W_%d[%d]: ", j, x);
            bn_print(msk->atts[j].w[x]);
        }
    }

}

void print_mpk(struct public_key_k_lin *mpk, const uint32_t n, const uint32_t ksec) {

    for (int i = 0; i < (ksec * (ksec + 1)); ++i) {
        printf("Matrix A_[%d]: \n", i);
        g1_print(mpk->a_mat[i]);
    }

    printf("\n");

    for (int x = 0; x < (n+1); ++x) {
        printf("Matrix AW_%x: \n", x);
        for (int y = 0; y < (ksec * ksec); ++y) {
            printf("Matrix AW_%d[%d]: \n", x, y);
            g1_print(mpk->mats[x].w[y]);
        }
        printf("\n");
    }

    printf("\n");

    for (int z = 0; z < ksec; ++z) {
        printf("Vector (mapping) e_[%d]: \n", z);
        gt_print(mpk->e_mat[z]);
    }


}

void print_ct(struct ciphertext_K_Lin *CT_A, const uint32_t n, const uint32_t ksec) {
    printf("Message M: \n");
    gt_print(CT_A->M);

    printf("\n");

    for (int i = 0; i < (ksec + 1); ++i) {
        printf("Vector ct_[%d]: \n", i);
        g1_print(CT_A->C_1[i]);
    }

    printf("\n");

    for (int x = 0; x < (n+1); ++x) {
        printf("vector ct_2,_%x: \n", x);
        for (int y = 0; y < ksec; ++y) {
            printf("Vector ct_2,_%d[%d]: \n", x, y);
            g1_print(CT_A->C_2[x].c_2_mat[y]);
        }
        printf("\n");
    }

    printf("Ct3: \n");
    gt_print(CT_A->C_3_one_val);
}

void test_matrix_mul_vector(const uint32_t ksec, bn_t o){
    bn_t one;
    init_null_new_bn_t_var(one);
    bn_set_dig(one, 55);

    bn_t five;
    init_null_new_bn_t_var(five);
    bn_set_dig(five, 130);

    bn_t eleven;
    init_null_new_bn_t_var(eleven);
    bn_set_dig(eleven, 110);

    printf("\n");


    const int matrixSize = (((2* ksec) + 1) * (ksec));
    const int vectorSize = ((2* ksec) + 1);

    bn_t matrix[matrixSize];
    bn_t vector[vectorSize];


    for (int i = 1; i < (matrixSize + 1); ++i) {
        init_null_new_bn_t_var(matrix[i-1]);
        bn_set_dig(matrix[i-1], i);
        bn_print(matrix[i-1]);
    }

    printf("\n");

    for (int i = 1; i < (vectorSize + 1); ++i) {
        init_null_new_bn_t_var(vector[i-1]);
        bn_set_dig(vector[i-1], i);
        bn_print(vector[i-1]);
    }

    printf("\n");

    bn_t *mat_vec;
    int a_rows = (ksec);
    int a_cols = ((2 * ksec) + 1);
    int v_rows = ((2 * ksec) + 1);

    bn_t out[ksec];

    mat_vec = matrix_mul_vector(out, matrix, vector, a_rows, a_cols, v_rows, o);

    for (int j = 0; j < a_rows; ++j) {
        printf("results \n");
        bn_print(mat_vec[j]);
    }

    printf("comparisons \n");
    bn_print(one);
    bn_print(five);
}

void test_vector_trans_mul_matrix(const uint32_t ksec, bn_t o) {
    bn_t two_two;
    init_null_new_bn_t_var(two_two);
    bn_set_dig(two_two, 7);

    bn_t two_eight;
    init_null_new_bn_t_var(two_eight);
    bn_set_dig(two_eight, 10);

    bn_t two_ten;
    init_null_new_bn_t_var(two_ten);
    bn_set_dig(two_ten, 110);

    printf("\n");

    const int vectorSize = (ksec);
    const int matrixSize = (ksec * (ksec));

    bn_t matrix[matrixSize];
    bn_t s_vector[vectorSize];

    for (int i = 1; i < (matrixSize + 1); ++i) {
        init_null_new_bn_t_var(matrix[i-1]);
        bn_set_dig(matrix[i-1], i);
        bn_print(matrix[i-1]);
    }

    printf("\n");

    for (int i = 1; i < (vectorSize + 1); ++i) {
        init_null_new_bn_t_var(s_vector[i-1]);
        bn_set_dig(s_vector[i-1], i);
        bn_print(s_vector[i-1]);
    }

    printf("\n");

    bn_t *ct_1;
    int a_cols = (ksec);
    int v_cols = (ksec);
    int a_rows = (ksec);

    bn_t out2[a_cols];
    ct_1 = vector_trans_mul_matrix(out2, s_vector, matrix, v_cols, a_cols, a_rows, o);

    for (int j = 0; j < a_cols; ++j) {
        printf("results \n");
        bn_print(ct_1[j]);
    }

    printf("comparisons \n");
    bn_print(two_two);
    bn_print(two_eight);
}

void test_matrix_mul_matrix(const uint32_t ksec, bn_t o) {
    bn_t two_two;
    init_null_new_bn_t_var(two_two);
    bn_set_dig(two_two, 38);

    bn_t two_eight;
    init_null_new_bn_t_var(two_eight);
    bn_set_dig(two_eight, 44);

    bn_t two_ten;
    init_null_new_bn_t_var(two_ten);
    bn_set_dig(two_ten, 50);

    bn_t two_twelve;
    init_null_new_bn_t_var(two_twelve);
    bn_set_dig(two_twelve, 56);

    bn_t six;
    init_null_new_bn_t_var(six);
    bn_set_dig(six, 83);

    bn_t fifteen;
    init_null_new_bn_t_var(fifteen);
    bn_set_dig(fifteen, 98);

    bn_t x1;
    init_null_new_bn_t_var(x1);
    bn_set_dig(x1, 113);

    bn_t x2;
    init_null_new_bn_t_var(x2);
    bn_set_dig(x2, 128);

    bn_t x3;
    init_null_new_bn_t_var(x3);
    bn_set_dig(x3, 128);

    bn_t x4;
    init_null_new_bn_t_var(x4);
    bn_set_dig(x4, 152);

    bn_t x5;
    init_null_new_bn_t_var(x5);
    bn_set_dig(x5, 176);

    bn_t x6;
    init_null_new_bn_t_var(x6);
    bn_set_dig(x6, 200);


    printf("\n");

    const int matrix1Size = ((ksec) * (ksec + 1));
    const int matrix2Size = ((ksec + 1) * (ksec));


    bn_t matrix1[matrix1Size];
    bn_t matrix2[matrix2Size];


    for (int i = 1; i < (matrix1Size + 1); ++i) {
        init_null_new_bn_t_var(matrix1[i-1]);
        bn_set_dig(matrix1[i-1], i);
        bn_print(matrix1[i-1]);
    }

    printf("\n");

    for (int i = 1; i < (matrix2Size + 1); ++i) {
        init_null_new_bn_t_var(matrix2[i-1]);
        bn_set_dig(matrix2[i-1], i);
        bn_print(matrix2[i-1]);
    }

    printf("\n");


    bn_t *ct_1;
    int m1_rows = (ksec);
    int m1_cols = (ksec + 1);
    int m2_rows = (ksec + 1);
    int m2_cols = (ksec);


    bn_t *m1_m2;
    bn_t out3[m1_rows * m2_cols];
    m1_m2 = matrixA_mul_matrixW(out3, matrix1, matrix2, m1_rows, m1_cols, m2_rows, m2_cols, o);

    for (int j = 0; j < (m1_rows * m2_cols); ++j) {
        printf("results \n");
        bn_print(m1_m2[j]);
    }

    printf("comparisons \n");
    bn_print(two_two);
    bn_print(two_eight);
    bn_print(two_ten);

    bn_print(two_twelve);
    bn_print(six);
    bn_print(fifteen);

    bn_print(x1);
    bn_print(x2);
    bn_print(x3);

    bn_print(x4);
    bn_print(x5);
    bn_print(x6);


    for (int w = 1; w <= 100; ++w) {
        bn_t tmp_dig;
        init_null_new_bn_t_var(tmp_dig);
        bn_set_dig(tmp_dig, w);
        printf("Bn_t value of the integer: %d: ", w);
        bn_print(tmp_dig);
    }

}

void test_vector_dot_product(const uint32_t ksec, bn_t o) {
    bn_t res;
    init_null_new_bn_t_var(res);
    bn_set_dig(res, 204);


    printf("\n");


    const int vector1size = (ksec + 6);
    const int vector2size = (ksec + 6);

    bn_t vector_1[vector1size];
    bn_t vector_2[vector2size];


    for (int i = 1; i < (vector1size + 1); ++i) {
        init_null_new_bn_t_var(vector_1[i-1]);
        bn_set_dig(vector_1[i-1], i);
        bn_print(vector_1[i-1]);
    }

    printf("\n");

    for (int i = 1; i < (vector2size + 1); ++i) {
        init_null_new_bn_t_var(vector_2[i-1]);
        bn_set_dig(vector_2[i-1], i);
        bn_print(vector_2[i-1]);
    }

    printf("\n");


    bn_t dot;
    init_null_new_bn_t_var(dot);
    vector_dot_product(dot, vector_1, vector_2, vector1size, vector2size, o);

    printf("results \n");
    bn_print(dot);

    printf("comparisons \n");
    bn_print(res);
}

void test_matrix_mul_scalar(const uint32_t ksec, bn_t o) {
    bn_t x1;
    init_null_new_bn_t_var(x1);
    bn_set_dig(x1, 2);

    bn_t x2;
    init_null_new_bn_t_var(x2);
    bn_set_dig(x2, 4);

    bn_t x3;
    init_null_new_bn_t_var(x3);
    bn_set_dig(x3, 6);

    bn_t x4;
    init_null_new_bn_t_var(x4);
    bn_set_dig(x4, 8);

    bn_t x5;
    init_null_new_bn_t_var(x5);
    bn_set_dig(x5, 10);

    bn_t x6;
    init_null_new_bn_t_var(x6);
    bn_set_dig(x6, 12);

    bn_t x7;
    init_null_new_bn_t_var(x7);
    bn_set_dig(x7, 14);

    bn_t x8;
    init_null_new_bn_t_var(x8);
    bn_set_dig(x8, 16);

    bn_t x9;
    init_null_new_bn_t_var(x9);
    bn_set_dig(x9, 18);

    bn_t x10;
    init_null_new_bn_t_var(x10);
    bn_set_dig(x10, 20);


    const int matrixSize = ((ksec) * ((2*ksec) + 1));
    bn_t matrix1[matrixSize];

    for (int i = 1; i < (matrixSize + 1); ++i) {
        init_null_new_bn_t_var(matrix1[i-1]);
        bn_set_dig(matrix1[i-1], i);
        bn_print(matrix1[i-1]);
    }

    printf("\n");

    int scalar = 2;
    printf("Scalar = %d \n", scalar);


    int m1_rows = ((2*ksec)+1);
    int m1_cols = (ksec);


    bn_t *m1_s;
    bn_t out3[matrixSize];
    m1_s = matrix_mul_scalar(out3, matrix1, scalar, m1_rows, m1_cols, o);

    for (int j = 0; j < (matrixSize); ++j) {
        printf("results \n");
        bn_print(m1_s[j]);
    }

    printf("comparisons \n");
    bn_print(x1);
    bn_print(x2);
    bn_print(x3);
    bn_print(x4);
    bn_print(x5);
    bn_print(x6);
    bn_print(x7);
    bn_print(x8);
    bn_print(x9);
    bn_print(x10);

}

void test_matrix_add_matrix(const uint32_t ksec, bn_t o) {
    bn_t x1;
    init_null_new_bn_t_var(x1);
    bn_set_dig(x1, 2);

    bn_t x2;
    init_null_new_bn_t_var(x2);
    bn_set_dig(x2, 4);

    bn_t x3;
    init_null_new_bn_t_var(x3);
    bn_set_dig(x3, 6);

    bn_t x4;
    init_null_new_bn_t_var(x4);
    bn_set_dig(x4, 8);

    bn_t x5;
    init_null_new_bn_t_var(x5);
    bn_set_dig(x5, 10);

    bn_t x6;
    init_null_new_bn_t_var(x6);
    bn_set_dig(x6, 12);

    bn_t x7;
    init_null_new_bn_t_var(x7);
    bn_set_dig(x7, 14);

    bn_t x8;
    init_null_new_bn_t_var(x8);
    bn_set_dig(x8, 16);

    bn_t x9;
    init_null_new_bn_t_var(x9);
    bn_set_dig(x9, 18);

    bn_t x10;
    init_null_new_bn_t_var(x10);
    bn_set_dig(x10, 20);


    const int matrix1Size = ((ksec) * ((2*ksec) + 1));
    const int matrix2Size = ((ksec) * ((2*ksec) + 1));
    bn_t matrix1[matrix1Size];
    bn_t matrix2[matrix2Size];

    for (int i = 1; i < (matrix1Size + 1); ++i) {
        init_null_new_bn_t_var(matrix1[i-1]);
        bn_set_dig(matrix1[i-1], i);
        bn_print(matrix1[i-1]);
    }

    printf("\n");

    for (int i = 1; i < (matrix2Size + 1); ++i) {
        init_null_new_bn_t_var(matrix2[i-1]);
        bn_set_dig(matrix2[i-1], i);
        bn_print(matrix2[i-1]);
    }


    int m1_rows = ((2*ksec)+1);
    int m1_cols = (ksec);
    int m2_rows = ((2*ksec)+1);
    int m2_cols = (ksec);


    bn_t *m1_add_m2;
    bn_t out3[matrix1Size];
    m1_add_m2 = matrix_add_matrix(out3, matrix1, matrix2, m1_rows, m1_cols, m2_rows, m2_cols, o);

    for (int j = 0; j < (matrix1Size); ++j) {
        printf("results \n");
        bn_print(m1_add_m2[j]);
    }

    printf("comparisons \n");
    bn_print(x1);
    bn_print(x2);
    bn_print(x3);
    bn_print(x4);
    bn_print(x5);
    bn_print(x6);
    bn_print(x7);
    bn_print(x8);
    bn_print(x9);
    bn_print(x10);
}

void test_vector_add_vector(const uint32_t ksec, bn_t o) {
    bn_t x1;
    init_null_new_bn_t_var(x1);
    bn_set_dig(x1, 2);

    bn_t x2;
    init_null_new_bn_t_var(x2);
    bn_set_dig(x2, 4);

    bn_t x3;
    init_null_new_bn_t_var(x3);
    bn_set_dig(x3, 6);

    bn_t x4;
    init_null_new_bn_t_var(x4);
    bn_set_dig(x4, 8);

    bn_t x5;
    init_null_new_bn_t_var(x5);
    bn_set_dig(x5, 10);


    const int vector1Size = (((2*ksec) + 1));
    const int vector2Size = (((2*ksec) + 1));
    bn_t vector1[vector1Size];
    bn_t vector2[vector2Size];

    for (int i = 1; i < (vector1Size + 1); ++i) {
        init_null_new_bn_t_var(vector1[i-1]);
        bn_set_dig(vector1[i-1], i);
        bn_print(vector1[i-1]);
    }

    printf("\n");

    for (int i = 1; i < (vector2Size + 1); ++i) {
        init_null_new_bn_t_var(vector2[i-1]);
        bn_set_dig(vector2[i-1], i);
        bn_print(vector2[i-1]);
    }


    bn_t *v1_add_v2;
    bn_t out3[vector1Size];
    v1_add_v2 = vector_add_vector(out3, vector1, vector2, vector1Size, vector2Size, o);

    for (int j = 0; j < (vector1Size); ++j) {
        printf("results \n");
        bn_print(v1_add_v2[j]);
    }

    printf("comparisons \n");
    bn_print(x1);
    bn_print(x2);
    bn_print(x3);
    bn_print(x4);
    bn_print(x5);
}

