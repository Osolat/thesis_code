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

#ifndef __STRUCTURES_H__
#define __STRUCTURES_H__
const int kss = 2;                                                    //Here kss=1 means (SXLIN) and kss=2 means (DLIN)

#include <stdint.h>

extern "C" {
#include <relic/relic.h>
}

// Key Policy

struct public_key_kp_gpsw {
    uint32_t N_ATTR;
    gt_t Y;
    g2_t *T_values;
};

struct public_key_kp_gpsw_oe {
    uint32_t N_ATTR;
    gt_t Y;
    g1_t *T_values;
};

struct master_key_kp_gpsw {
    uint32_t N_ATTR;
    struct attribute *attributes;
    bn_t y;
    bn_t *t_values;
};

struct secret_key_kp_gpsw {
    g1_t *D_values;
};

struct secret_key_kp_gpsw_oe {
    g2_t *D_values;
};

struct ciphertext_kp_gpsw {
    bn_t gamma;
    gt_t E_prime;
    g2_t *E_values;
};

struct ciphertext_kp_gpsw_oe {
    bn_t gamma;
    gt_t E_prime;
    g1_t *E_values;
};

int init_master_key_kp_gpsw(const uint32_t n_attr, struct master_key_kp_gpsw *m);
int init_public_key_kp_gpsw(const uint32_t n_attr, struct public_key_kp_gpsw *p);
int init_public_key_kp_gpsw_oe(const uint32_t n_attr, struct public_key_kp_gpsw_oe *p);
int init_secret_key_kp_gpsw(const uint32_t n_attr, struct secret_key_kp_gpsw *sk);
int init_secret_key_kp_gpsw_oe(const uint32_t n_attr, struct secret_key_kp_gpsw_oe *sk);

int init_ciphertext_kp_gpsw(const uint32_t n_attr, struct ciphertext_kp_gpsw *E);
int init_ciphertext_kp_gpsw_oe(const uint32_t n_attr, struct ciphertext_kp_gpsw_oe *E);

struct master_key_kp_gpsw_lu {
    bn_t y;
};

struct public_key_kp_gpsw_lu {
    g1_t g1;
    g2_t g2;
    // TODO: Most likely don't need t_values if using T as a hashing function -> G_2
    g2_t *t_values;
};

struct secret_key_kp_gpsw_lu {
    g2_t *D_values;
    g1_t *R_values;
};

struct ciphertext_kp_gpsw_lu {
    bn_t gamma;
    gt_t E_prime;
    g1_t E_prime_prime;
    g2_t *E_values;
};

int init_master_key_kp_gpsw_lu(const uint32_t n_attr, struct master_key_kp_gpsw_lu *m);
int init_public_key_kp_gpsw_lu(const uint32_t n_attr, struct public_key_kp_gpsw_lu *p);
int init_secret_key_kp_gpsw_lu(const uint32_t n_attr, struct secret_key_kp_gpsw_lu *sk);
int init_ciphertext_kp_gpsw_lu(const uint32_t n_attr, struct ciphertext_kp_gpsw_lu *E);


//Stuff used for the K_LIN scheme.
//const int k = 2;                                     //Here k=1 means (SXLIN) and k=2 means (DLIN)

//TODO factor out the constant k such that it is set at compile time.
struct k_lin_att {
    uint32_t attr;
    bn_t w[(kss+1)*kss];                                 //Matrix w has size ((k+1) * k)
};

struct k_lin_mat {
    uint32_t attr;
    g1_t w[kss*kss];                                      //Matrix AW_i will have size k*k after matrix multiplication.
};

struct k_lin_secret_key {
    uint32_t attr;
    g2_t sk_one[kss+1];                                //SK1 will be a vector of k+1 since matrix-vector multiplication yields a vector of size k. Then addition of two vectors must have same dimensions, so we finally get a vector of size k+1
    g2_t sk_two[kss];                                  //SK2 is a vector of size k.
};

struct tmp_vj {
    bn_t vec_j[kss+1];                                 //The vj-vector will have size k+1 otherwise we can't do vector addition for SK1.
};

struct tmp_rj {
    bn_t vec_rj[kss];                                  //The rj-vector has size k.
};

struct c_attribute_K_Lin {
    uint32_t attr;
    g1_t c_2_mat[kss];                                 //Each CT_2 will be a vector of size k, since AWi yields a k*k matrix and that multiplied with a transposed vector of size k will yield a vector of size k.
};

struct master_key_k_lin {
    uint32_t N_ATTR;
    uint32_t N_SEC;
    struct k_lin_att *atts;                          //Set of matrices W_i, where the set-size is N_ATTR + 1, since W_0 = 0 and W_1,...,W_n is valid matrices of msk.
    bn_t *v_share;                                   //A vector v of size (k+1).
};

struct public_key_k_lin {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    struct k_lin_mat *mats;                          //Set of matrices AW_i, where the set-size is N_ATTR + 1, since AW_0 = 0 by matrix multiplication of W_0 and AW_1,...,AW_n is valid matrices of mpk.
    g1_t *a_mat;                                     //The matrix A has size k*(k+1) and the entries are group1 elements raised to the power of the matrix-product AW_i of each component.
    gt_t *e_mat;                                     //Asymmetric mapping E = e([A]_1, [v]_2) = gt^Av = e(g,h)^Av which has size k, since Av yields a vector of size k. //TODO can't malloc gt_t as it gives a bus core seg fault.
};

struct secret_key_K_Lin {
    uint32_t N_ATTR;
    struct k_lin_secret_key * sk;                    //Set of secret key pairs (SK1,SK2) where sk1 is a vector of size k+1 and SK2 is a vector of size k.
};

struct sk_tmp_vj {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    struct tmp_vj * vj;                              //Temporary set of vj vectors containing the j'th share of vector v. The set-size is j, where j=N_ATTR when the policy-tree only consists of AND gates. //TODO adjust for policy-trees with or gates as well.
    struct tmp_rj * rj;                              //Temporary set of rj vectors containing random bn_t values. The set-size is j, where j=N_ATTR when the policy-tree only consists of AND gates. //TODO adjust for policy-trees with or gates as well.
};

struct ciphertext_K_Lin {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    gt_t M;                                         //This is the message M from the encryption algorithm, sets it to 1 for simplicity.
    g1_t *C_1;                                      //CT1 is a vector of size k+1 because the vector-matrix multiplication s^T*A yields a vector of size k+1
    struct c_attribute_K_Lin * C_2;                 //CT2 is a set of vectors where the set-size is based on #attributes where x_i = 1?
    gt_t C_3_one_val;                               //CT3 is a vector of size k+1 because of the cross-product of two vectors of size k+1 yields a vector of size k+1 //TODO cannot change like for C1 since it give a bus seg fault

};

int init_ciphertext_K_Lin(const uint32_t n_attr, const uint32_t kss, struct ciphertext_K_Lin *c);
int init_secret_key_K_Lin(const uint32_t n_attr, struct secret_key_K_Lin *s);
int init_sk_tmp_vj(const uint32_t n_attr, const uint32_t kss, struct sk_tmp_vj *v);
int init_master_key_k_lin(const uint32_t n_attr, const uint32_t kss, struct master_key_k_lin *m);
int init_public_key_k_lin(const uint32_t n_attr, const uint32_t kss, struct public_key_k_lin *p);


//Stuff used for KLin large universe



struct c_attribute_K_Lin_lu_c23 {
    uint32_t attr;
    g1_t c_2_vec[kss];                                                            //Each CT_2 will be a vector of size k, since AWi yields a k*k matrix and that multiplied with a transposed vector of size k will yield a vector of size k.
    g1_t c_3_vec[(2*kss)+1];
};


struct tmp_vj_lu {
    bn_t vec_j[(2*kss)+1];                                                           //The vj-vector will have size 2k+1 otherwise we can't do vector addition for SK1.
};

struct tmp_rj_lu {
    bn_t vec_rj[kss];                                                             //The rj-vector has size k.
};

struct si_lu {
    bn_t si_vec[kss];                                                             //The rj-vector has size k.
};

struct sk_tmp_vectors_lu {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    struct tmp_vj_lu * vj;                                                     //Temporary set of vj vectors containing the j'th share of vector v. The set-size is j, where j=N_ATTR when the policy-tree only consists of AND gates. //TODO adjust for policy-trees with or gates as well.
    struct tmp_rj_lu * rj;                                                     //Temporary set of rj vectors containing random bn_t values. The set-size is j, where j=N_ATTR when the policy-tree only consists of AND gates. //TODO adjust for policy-trees with or gates as well.
};

struct tmp_si_lu {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    struct si_lu * si;                                                     //Temporary set of vj vectors containing the j'th share of vector v. The set-size is j, where j=N_ATTR when the policy-tree only consists of AND gates. //TODO adjust for policy-trees with or gates as well.
};



struct k_lin_secret_key_lu_13 {
    uint32_t attr;
    g2_t sk_one[(2*kss)+1];                                                      //sk1 is a vector of size 2k+1
    g2_t sk_two[kss];                                                         //sk2 is a vector of size k
    g2_t sk_three[(2*kss)+1];                                                    //sk3 is a vector of size 2k+1 not exactly sure what rho(j)!=0 is and how it works for pure and gates
};

struct k_lin_secret_key_lu_4 {
    uint32_t attr;
    g2_t sk_four[(2*kss)+1];                                                     //sk4 is a vector of size 2k+1 and again not sure about rho(j)=0.
};



struct master_key_k_lin_lu {
    uint32_t N_ATTR;
    uint32_t N_SEC;
    bn_t *W_matrix;                                                          //w matrix of size (2k+1)*k
    bn_t *W0_matrix;                                                         //w0 matrix of size (2k+1)*k
    bn_t *W1_matrix;                                                         //w1 matrix of size (2k+1)*k
    bn_t *v_secret;                                                          //vector of secrets  of size 2k+1
};


struct public_key_k_lin_lu {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    g1_t *A1_mat;                                                           //A1 matrix of size K*(2k+1)
    g1_t *AW_mat;                                                           //AW matrix of size k*K
    g1_t *AW0_mat;                                                          //AW0 matrix of size k*K
    g1_t *AW1_mat;                                                          //AW1 matrix of size k*K
    gt_t *e_mat;                                                            //e mapping vector of size k.
};


struct secret_key_K_Lin_lu {
    uint32_t N_ATTR;
    struct k_lin_secret_key_lu_13 * sk13;
    struct k_lin_secret_key_lu_4 * sk4;
};


struct ciphertext_K_Lin_lu {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    gt_t M;                                                             //This is the message M from the encryption algorithm, sets it to 1 for simplicity.
    g1_t *C_1;                                                          //CT1 is a vector of size 2k+1 because the vector-matrix multiplication s^T*A yields a vector of size 2k+1
    struct c_attribute_K_Lin_lu_c23 * C_23;                            //CT2 is a set of vectors where the set-size is based on #attributes where x_i = 1?
    gt_t C_4_one_val;                                                   //CT4 one gt value.

};

int init_ciphertext_K_Lin_lu(const uint32_t n_attr, const uint32_t kss, struct ciphertext_K_Lin_lu *c);
int init_secret_key_K_Lin_lu(const uint32_t n_attr, struct secret_key_K_Lin_lu *s);
int init_sk_tmp_vectors_lu(const uint32_t n_attr, const uint32_t kss, struct sk_tmp_vectors_lu *v);
int init_master_key_k_lin_lu(const uint32_t n_attr, const uint32_t kss, struct master_key_k_lin_lu *m);
int init_public_key_k_lin_lu(const uint32_t n_attr, const uint32_t kss, struct public_key_k_lin_lu *p);
int init_tmp_si_lu(const uint32_t n_attr, const uint32_t kss, struct tmp_si_lu *si);


//K-Lin unbounded OK structures
struct c_attribute_K_Lin_ub_c23_ok {
    uint32_t attr;
    g2_t c_2_vec[kss];                                                            //Each CT_2 will be a vector of size k, since AWi yields a k*k matrix and that multiplied with a transposed vector of size k will yield a vector of size k.
    g2_t c_3_vec[(2*kss)+1];
};


struct tmp_vj_ub_ok {
    bn_t vec_j[(2*kss)+1];                                                           //The vj-vector will have size 2k+1 otherwise we can't do vector addition for SK1.
};

struct tmp_rj_ub_ok {
    bn_t vec_rj[kss];                                                             //The rj-vector has size k.
};

struct si_ub_ok {
    bn_t si_vec[kss];                                                             //The rj-vector has size k.
};

struct sk_tmp_vectors_ub_ok {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    struct tmp_vj_ub_ok * vj;                                                     //Temporary set of vj vectors containing the j'th share of vector v. The set-size is j, where j=N_ATTR when the policy-tree only consists of AND gates. //TODO adjust for policy-trees with or gates as well.
    struct tmp_rj_ub_ok * rj;                                                     //Temporary set of rj vectors containing random bn_t values. The set-size is j, where j=N_ATTR when the policy-tree only consists of AND gates. //TODO adjust for policy-trees with or gates as well.
};

struct tmp_si_ub_ok {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    struct si_ub_ok * si;                                                     //Temporary set of vj vectors containing the j'th share of vector v. The set-size is j, where j=N_ATTR when the policy-tree only consists of AND gates. //TODO adjust for policy-trees with or gates as well.
};



struct k_lin_secret_key_ub_13_ok {
    uint32_t attr;
    g1_t sk_one[(2*kss)+1];                                                      //sk1 is a vector of size 2k+1
    g1_t sk_two[kss];                                                         //sk2 is a vector of size k
    g1_t sk_three[(2*kss)+1];                                                    //sk3 is a vector of size 2k+1 not exactly sure what rho(j)!=0 is and how it works for pure and gates
};

struct k_lin_secret_key_ub_4_ok {
    uint32_t attr;
    g1_t sk_four[(2*kss)+1];                                                     //sk4 is a vector of size 2k+1 and again not sure about rho(j)=0.
};



struct master_key_k_lin_ub_ok {
    uint32_t N_ATTR;
    uint32_t N_SEC;
    bn_t *W_matrix;                                                          //w matrix of size (2k+1)*k
    bn_t *W0_matrix;                                                         //w0 matrix of size (2k+1)*k
    bn_t *W1_matrix;                                                         //w1 matrix of size (2k+1)*k
    bn_t *v_secret;                                                          //vector of secrets  of size 2k+1
};


struct public_key_k_lin_ub_ok {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    g2_t *A1_mat;                                                           //A1 matrix of size K*(2k+1)
    g2_t *AW_mat;                                                           //AW matrix of size k*K
    g2_t *AW0_mat;                                                          //AW0 matrix of size k*K
    g2_t *AW1_mat;                                                          //AW1 matrix of size k*K
    gt_t *e_mat;                                                            //e mapping vector of size k.
};


struct secret_key_K_Lin_ub_ok {
    uint32_t N_ATTR;
    struct k_lin_secret_key_ub_13_ok * sk13;
    struct k_lin_secret_key_ub_4_ok * sk4;
};


struct ciphertext_K_Lin_ub_ok {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    gt_t M;                                                             //This is the message M from the encryption algorithm, sets it to 1 for simplicity.
    g2_t *C_1;                                                          //CT1 is a vector of size 2k+1 because the vector-matrix multiplication s^T*A yields a vector of size 2k+1
    struct c_attribute_K_Lin_ub_c23_ok * C_23;                            //CT2 is a set of vectors where the set-size is based on #attributes where x_i = 1?
    gt_t C_4_one_val;                                                   //CT4 one gt value.

};

int init_ciphertext_K_Lin_ub_ok(const uint32_t n_attr, const uint32_t kss, struct ciphertext_K_Lin_ub_ok *c);
int init_secret_key_K_Lin_ub_ok(const uint32_t n_attr, struct secret_key_K_Lin_ub_ok *s);
int init_sk_tmp_vectors_ub_ok(const uint32_t n_attr, const uint32_t kss, struct sk_tmp_vectors_ub_ok *v);
int init_master_key_k_lin_ub_ok(const uint32_t n_attr, const uint32_t kss, struct master_key_k_lin_ub_ok *m);
int init_public_key_k_lin_ub_ok(const uint32_t n_attr, const uint32_t kss, struct public_key_k_lin_ub_ok *p);
int init_tmp_si_ub_ok(const uint32_t n_attr, const uint32_t kss, struct tmp_si_ub_ok *si);


//K-Lin unbounded OD structures
struct c_attribute_K_Lin_ub_c23_od {
    uint32_t attr;
    g1_t c_2_vec[kss];                                                            //Each CT_2 will be a vector of size k, since AWi yields a k*k matrix and that multiplied with a transposed vector of size k will yield a vector of size k.
    g1_t c_3_vec[(2*kss)+1];
};


struct tmp_vj_ub_od {
    bn_t vec_j[(2*kss)+1];                                                           //The vj-vector will have size 2k+1 otherwise we can't do vector addition for SK1.
};

struct tmp_rj_ub_od {
    bn_t vec_rj[kss];                                                             //The rj-vector has size k.
};

struct si_ub_od {
    bn_t si_vec[kss];                                                             //The rj-vector has size k.
};

struct sk_tmp_vectors_ub_od {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    struct tmp_vj_ub_od * vj;                                                     //Temporary set of vj vectors containing the j'th share of vector v. The set-size is j, where j=N_ATTR when the policy-tree only consists of AND gates. //TODO adjust for policy-trees with or gates as well.
    struct tmp_rj_ub_od * rj;                                                     //Temporary set of rj vectors containing random bn_t values. The set-size is j, where j=N_ATTR when the policy-tree only consists of AND gates. //TODO adjust for policy-trees with or gates as well.
};

struct tmp_si_ub_od {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    struct si_ub_od * si;                                                     //Temporary set of vj vectors containing the j'th share of vector v. The set-size is j, where j=N_ATTR when the policy-tree only consists of AND gates. //TODO adjust for policy-trees with or gates as well.
};



struct k_lin_secret_key_ub_13_od {
    uint32_t attr;
    g1_t sk_one[(2*kss)+1];                                                      //sk1 is a vector of size 2k+1
    g2_t sk_two[kss];                                                         //sk2 is a vector of size k
    g2_t sk_three[(2*kss)+1];                                                    //sk3 is a vector of size 2k+1 not exactly sure what rho(j)!=0 is and how it works for pure and gates
};

struct k_lin_secret_key_ub_4_od {
    uint32_t attr;
    g1_t sk_four[(2*kss)+1];                                                     //sk4 is a vector of size 2k+1 and again not sure about rho(j)=0.
};



struct master_key_k_lin_ub_od {
    uint32_t N_ATTR;
    uint32_t N_SEC;
    bn_t *W_matrix;                                                          //w matrix of size (2k+1)*k
    bn_t *W0_matrix;                                                         //w0 matrix of size (2k+1)*k
    bn_t *W1_matrix;                                                         //w1 matrix of size (2k+1)*k
    bn_t *v_secret;                                                          //vector of secrets  of size 2k+1
};


struct public_key_k_lin_ub_od {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    g1_t *A1_mat_g1;
    g2_t *A1_mat_g2;                                                           //A1 matrix of size K*(2k+1)
    g1_t *AW_mat;                                                           //AW matrix of size k*K
    g1_t *AW0_mat;                                                          //AW0 matrix of size k*K
    g1_t *AW1_mat;                                                          //AW1 matrix of size k*K
    gt_t *e_mat;                                                            //e mapping vector of size k.
};


struct secret_key_K_Lin_ub_od {
    uint32_t N_ATTR;
    struct k_lin_secret_key_ub_13_od * sk13;
    struct k_lin_secret_key_ub_4_od * sk4;
};


struct ciphertext_K_Lin_ub_od {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    gt_t M;                                                             //This is the message M from the encryption algorithm, sets it to 1 for simplicity.
    g2_t *C_1;                                                          //CT1 is a vector of size 2k+1 because the vector-matrix multiplication s^T*A yields a vector of size 2k+1
    struct c_attribute_K_Lin_ub_c23_od * C_23;                            //CT2 is a set of vectors where the set-size is based on #attributes where x_i = 1?
    gt_t C_4_one_val;                                                   //CT4 one gt value.

};

int init_ciphertext_K_Lin_ub_od(const uint32_t n_attr, const uint32_t kss, struct ciphertext_K_Lin_ub_od *c);
int init_secret_key_K_Lin_ub_od(const uint32_t n_attr, struct secret_key_K_Lin_ub_od *s);
int init_sk_tmp_vectors_ub_od(const uint32_t n_attr, const uint32_t kss, struct sk_tmp_vectors_ub_od *v);
int init_master_key_k_lin_ub_od(const uint32_t n_attr, const uint32_t kss, struct master_key_k_lin_ub_od *m);
int init_public_key_k_lin_ub_od(const uint32_t n_attr, const uint32_t kss, struct public_key_k_lin_ub_od *p);
int init_tmp_si_ub_od(const uint32_t n_attr, const uint32_t kss, struct tmp_si_ub_od *si);


struct k_lin_att_ok {
    uint32_t attr;
    bn_t w[(kss+1)*kss];                                 //Matrix w has size ((k+1) * k)
};

struct k_lin_mat_ok {
    uint32_t attr;
    g2_t w[kss*kss];                                      //Matrix AW_i will have size k*k after matrix multiplication.
};

struct k_lin_secret_key_ok {
    uint32_t attr;
    g1_t sk_one[kss+1];                                //SK1 will be a vector of k+1 since matrix-vector multiplication yields a vector of size k. Then addition of two vectors must have same dimensions, so we finally get a vector of size k+1
    g1_t sk_two[kss];                                  //SK2 is a vector of size k.
};

struct tmp_vj_ok {
    bn_t vec_j[kss+1];                                 //The vj-vector will have size k+1 otherwise we can't do vector addition for SK1.
};

struct tmp_rj_ok {
    bn_t vec_rj[kss];                                  //The rj-vector has size k.
};

struct c_attribute_K_Lin_ok {
    uint32_t attr;
    g2_t c_2_mat[kss];                                 //Each CT_2 will be a vector of size k, since AWi yields a k*k matrix and that multiplied with a transposed vector of size k will yield a vector of size k.
};

struct master_key_k_lin_ok {
    uint32_t N_ATTR;
    uint32_t N_SEC;
    struct k_lin_att_ok *atts;                          //Set of matrices W_i, where the set-size is N_ATTR + 1, since W_0 = 0 and W_1,...,W_n is valid matrices of msk.
    bn_t *v_share;                                   //A vector v of size (k+1).
};

struct public_key_k_lin_ok {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    struct k_lin_mat_ok *mats;                          //Set of matrices AW_i, where the set-size is N_ATTR + 1, since AW_0 = 0 by matrix multiplication of W_0 and AW_1,...,AW_n is valid matrices of mpk.
    g2_t *a_mat;                                     //The matrix A has size k*(k+1) and the entries are group1 elements raised to the power of the matrix-product AW_i of each component.
    gt_t *e_mat;                                     //Asymmetric mapping E = e([A]_1, [v]_2) = gt^Av = e(g,h)^Av which has size k, since Av yields a vector of size k. //TODO can't malloc gt_t as it gives a bus core seg fault.
};

struct secret_key_K_Lin_ok {
    uint32_t N_ATTR;
    struct k_lin_secret_key_ok * sk;                    //Set of secret key pairs (SK1,SK2) where sk1 is a vector of size k+1 and SK2 is a vector of size k.
};

struct sk_tmp_vj_ok {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    struct tmp_vj_ok* vj;                              //Temporary set of vj vectors containing the j'th share of vector v. The set-size is j, where j=N_ATTR when the policy-tree only consists of AND gates. //TODO adjust for policy-trees with or gates as well.
    struct tmp_rj_ok* rj;                              //Temporary set of rj vectors containing random bn_t values. The set-size is j, where j=N_ATTR when the policy-tree only consists of AND gates. //TODO adjust for policy-trees with or gates as well.
};

struct ciphertext_K_Lin_ok {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    gt_t M;                                         //This is the message M from the encryption algorithm, sets it to 1 for simplicity.
    g2_t *C_1;                                      //CT1 is a vector of size k+1 because the vector-matrix multiplication s^T*A yields a vector of size k+1
    struct c_attribute_K_Lin_ok * C_2;                 //CT2 is a set of vectors where the set-size is based on #attributes where x_i = 1?
    gt_t C_3_one_val;                               //CT3 is a vector of size k+1 because of the cross-product of two vectors of size k+1 yields a vector of size k+1 //TODO cannot change like for C1 since it give a bus seg fault

};

int init_ciphertext_K_Lin_ok(const uint32_t n_attr, const uint32_t kss, struct ciphertext_K_Lin_ok *c);
int init_secret_key_K_Lin_ok(const uint32_t n_attr, struct secret_key_K_Lin_ok *s);
int init_sk_tmp_vj_ok(const uint32_t n_attr, const uint32_t kss, struct sk_tmp_vj_ok *v);
int init_master_key_k_lin_ok(const uint32_t n_attr, const uint32_t kss, struct master_key_k_lin_ok *m);
int init_public_key_k_lin_ok(const uint32_t n_attr, const uint32_t kss, struct public_key_k_lin_ok *p);


//Stuff for k-lin_od
//TODO move stuff in right groups, most importantly move ct1 to h and sk1 to g. May require other changes of group.
struct k_lin_att_od {
    uint32_t attr;
    bn_t w[(kss+1)*kss];                                 //Matrix w has size ((k+1) * k)
};

struct k_lin_mat_od {
    uint32_t attr;
    g1_t w[kss*kss];                                      //Matrix AW_i will have size k*k after matrix multiplication.
};

struct k_lin_secret_key_od {
    uint32_t attr;
    g1_t sk_one[kss+1];                                //SK1 will be a vector of k+1 since matrix-vector multiplication yields a vector of size k. Then addition of two vectors must have same dimensions, so we finally get a vector of size k+1
    g2_t sk_two[kss];                                  //SK2 is a vector of size k.
};

struct tmp_vj_od {
    bn_t vec_j[kss+1];                                 //The vj-vector will have size k+1 otherwise we can't do vector addition for SK1.
};

struct tmp_rj_od {
    bn_t vec_rj[kss];                                  //The rj-vector has size k.
};

struct c_attribute_K_Lin_od {
    uint32_t attr;
    g1_t c_2_mat[kss];                                 //Each CT_2 will be a vector of size k, since AWi yields a k*k matrix and that multiplied with a transposed vector of size k will yield a vector of size k.
};

struct master_key_k_lin_od {
    uint32_t N_ATTR;
    uint32_t N_SEC;
    struct k_lin_att_od *atts;                          //Set of matrices W_i, where the set-size is N_ATTR + 1, since W_0 = 0 and W_1,...,W_n is valid matrices of msk.
    bn_t *v_share;                                   //A vector v of size (k+1).
};

struct public_key_k_lin_od {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    struct k_lin_mat_od *mats;                          //Set of matrices AW_i, where the set-size is N_ATTR + 1, since AW_0 = 0 by matrix multiplication of W_0 and AW_1,...,AW_n is valid matrices of mpk.
    g2_t *a_mat;                                     //The matrix A has size k*(k+1) and the entries are group1 elements raised to the power of the matrix-product AW_i of each component.
    gt_t *e_mat;                                     //Asymmetric mapping E = e([A]_1, [v]_2) = gt^Av = e(g,h)^Av which has size k, since Av yields a vector of size k. //TODO can't malloc gt_t as it gives a bus core seg fault.
};

struct secret_key_K_Lin_od {
    uint32_t N_ATTR;
    struct k_lin_secret_key_od * sk;                    //Set of secret key pairs (SK1,SK2) where sk1 is a vector of size k+1 and SK2 is a vector of size k.
};

struct sk_tmp_vj_od {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    struct tmp_vj_od * vj;                              //Temporary set of vj vectors containing the j'th share of vector v. The set-size is j, where j=N_ATTR when the policy-tree only consists of AND gates. //TODO adjust for policy-trees with or gates as well.
    struct tmp_rj_od * rj;                              //Temporary set of rj vectors containing random bn_t values. The set-size is j, where j=N_ATTR when the policy-tree only consists of AND gates. //TODO adjust for policy-trees with or gates as well.
};

struct ciphertext_K_Lin_od {
    uint32_t N_ATTR;
    uint32_t K_SEC;
    gt_t M;                                         //This is the message M from the encryption algorithm, sets it to 1 for simplicity.
    g2_t *C_1;                                      //CT1 is a vector of size k+1 because the vector-matrix multiplication s^T*A yields a vector of size k+1
    struct c_attribute_K_Lin_od * C_2;                 //CT2 is a set of vectors where the set-size is based on #attributes where x_i = 1?
    gt_t C_3_one_val;                               //CT3 is a vector of size k+1 because of the cross-product of two vectors of size k+1 yields a vector of size k+1 //TODO cannot change like for C1 since it give a bus seg fault

};

int init_ciphertext_K_Lin_od(const uint32_t n_attr, const uint32_t kss, struct ciphertext_K_Lin_od *c);
int init_secret_key_K_Lin_od(const uint32_t n_attr, struct secret_key_K_Lin_od *s);
int init_sk_tmp_vj_od(const uint32_t n_attr, const uint32_t kss, struct sk_tmp_vj_od *v);
int init_master_key_k_lin_od(const uint32_t n_attr, const uint32_t kss, struct master_key_k_lin_od *m);
int init_public_key_k_lin_od(const uint32_t n_attr, const uint32_t kss, struct public_key_k_lin_od *p);



#endif
