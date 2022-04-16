#include <cstdio>
#include <string>
#include <cmath>
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

int powint(int base, int exp)
{
	assert(exp > 0 && "powint: exp parameter has negative value");

	int result{ 1 };
	while (exp)
	{
		if (exp & 1)
			result *= base;
		exp >>= 1;
		base *= base;
	}

	return result;
}

//Thanks to Titas Chanda from stack overflow 
int coeff(int k,const vector<int>& roots)
{

  int size = roots.size(); // total no. of roots

  int loop_no = k; // total no. of nested loops 
  vector<int> loop_counter(loop_no+1); // loop_counter's are the actual iterators through the nested loops
  // like i_1, i_2, i_3 etc., loop_counter[loop_no] is needed to maintain the condition of the loops

  for(int i=0; i<loop_no+1; ++i)
    loop_counter[i]=0;   // initialize all iterators to zero


  vector<int> MAX(loop_no+1); // MAX value for a loop, like  for(int i_1=0; i_1 < MAX[1]; i++)... 
    for(int i=0 ; i<loop_no ; ++i)
      MAX[i] = size - loop_no  + i + 1; // calculated from the algorithm

    MAX[loop_no]=2; // do'es no matter, just != 1

    int  p1=0; // required to maintain the condition of the loops

    int sum(0); // sum of all products

    while(loop_counter[loop_no]==0)
    {
      // variable nested loops starts here

      int counter(0);
      // check that i_1 < i_2 < i_3 ....
      for(int i = 1 ; i < loop_no; ++i)
      {
        if(loop_counter[i-1] < loop_counter[i])
          ++counter;
      }


      if(counter == loop_no - 1) // true if i_1 < i_2 < i_3 ....
      {
        int prod(1);
        for(int i = 0 ; i < loop_no ; ++i)  
          prod *= roots[loop_counter[i]];   // taking products

        sum += prod;  // increament
      }


    // variable nested loops ends here...


    ++loop_counter[0];
    while(loop_counter[p1]==MAX[p1])
      {
        loop_counter[p1]=0;
        loop_counter[++p1]++;
        if(loop_counter[p1]!=MAX[p1])
          p1=0;
      }
    }
    return powint(-1,k)*sum;   // DO NOT FORGET THE NEGATIVE SIGN

}
//Inductive coefficients from roots computation
void ind_coeff(vector<vector<int>>& coeffs, const vector<int> roots, int n) {
    for (auto i = 0; i < n+1; i++){
        coeffs[i][i] = 1;
    }
    for (int i = 1; i < n+1; i++){
        for (int j = 0; j < i; j++) {
            if (j == 0) {
                coeffs[i][j] = -1*coeffs[i-1][0]*roots[i-1];
            } else {
                coeffs[i][j] = -1*coeffs[i-1][j]*roots[i-1]+coeffs[i-1][j-1];
            }
        }
    }
}
void ind_coeff_bn(bn_t *coeffs, const bn_t *roots, int n) {
    for (int i = 0; i < n; i++) {
        bn_null(coeffs[n*i+i]); bn_new(coeffs[n*i+i]);
        bn_set_dig(coeffs[n*i+i], 1);
    }
    for (int i = 1; i < n; i++){
        for (int j = 0; j < i; j++) {
            bn_null(coeffs[n*i+j]); bn_new(coeffs[n*i+j]);
            if (j == 0) {  
                bn_t tmp; bn_null(tmp); bn_new(tmp);
                bn_mul(tmp, coeffs[n*(i-1)], roots[i-1]);
                bn_neg(coeffs[n*i+j], tmp);
            } else {
                bn_t tmp; bn_null(tmp); bn_new(tmp);
                bn_mul(tmp, coeffs[n*(i-1)+j], roots[i-1]);
                bn_neg(tmp, tmp);
                bn_add(coeffs[n*i+j], tmp, coeffs[n*(i-1)+j-1]);
            }
        }
    } 
}


unsigned long long t[NTESTS];

std::vector<policy_coefficient> lsss_vector;

void setup_naive_oe(struct alp_pp_naive_oe *pp, bn_t alpha, bn_t order, int bound) {
    g1_t g1; g1_null(g1); g1_new(g1);
    g2_t g2; g2_null(g2); g2_new(g2);

    g1_get_gen(g1);
    //cout << "g1\n";
    //g1_print(g1);
    g2_get_gen(g2);
    //cout << "g2\n";
    //g2_print(g2);
    init_public_params_alp_oe(bound, pp);     
    pp->bound = bound;
    
    for(int i = 0; i < bound+1; i++) {
        bn_t alpha_i; bn_null(alpha_i); bn_new(alpha_i);
        g1_null(pp->U1[i]); g1_new(pp->U1[i]);
        g2_null(pp->U2[i]); g1_new(pp->U2[i]);
        bn_rand_mod(alpha_i, order); 
        g1_mul(pp->U1[i], g1, alpha_i);
        g2_mul(pp->U2[i], g2, alpha_i);
        if (i < bound) {
            bn_t beta_i; bn_null(beta_i); bn_new(beta_i);
            g1_null(pp->H1[i]); g1_new(pp->H1[i]);
            g2_null(pp->H2[i]); g2_new(pp->H2[i]);
            bn_rand_mod(beta_i, order);
            g1_mul(pp->H1[i], g1, beta_i);
            g2_mul(pp->H2[i], g2, beta_i);
        }
    }

    g1_copy(pp->g1, g1);
    g2_copy(pp->g2, g2);
    bn_copy(pp->order, order);
    pc_map(pp->gt, g1, g2);
    gt_exp(pp->gt, pp->gt, alpha);
}

void setup_g_oe(struct alp_pp_naive_oe *pp, bn_t alpha, bn_t order, int bound) {
    g1_t g1; g1_null(g1); g1_new(g1);
    g2_t g2; g2_null(g2); g2_new(g2);

    g1_get_gen(g1);
    //cout << "g1\n";
    //g1_print(g1);
    g2_get_gen(g2);
    //cout << "g2\n";
    //g2_print(g2);
    init_public_params_alp_oe(bound, pp);     
    pp->bound = bound;
    
    for(int i = 0; i < bound+1; i++) {
        bn_t alpha_i; bn_null(alpha_i); bn_new(alpha_i);
        g1_null(pp->U1[i]); g1_new(pp->U1[i]);
        g2_null(pp->U2[i]); g1_new(pp->U2[i]);
        bn_rand_mod(alpha_i, order); 
        g1_mul(pp->U1[i], g1, alpha_i);
        g2_mul(pp->U2[i], g2, alpha_i);
        if (i < bound) {
            bn_t beta_i; bn_null(beta_i); bn_new(beta_i);
            g1_null(pp->H1[i]); g1_new(pp->H1[i]);
            g2_null(pp->H2[i]); g2_new(pp->H2[i]);
            bn_rand_mod(beta_i, order);
            g1_mul(pp->H1[i], g1, beta_i);
            g2_mul(pp->H2[i], g2, beta_i);
        }
    }
    g1_copy(pp->g1, g1);
    g2_copy(pp->g2, g2);
    bn_copy(pp->order, order);
    g1_t g_alpha; g1_null(g_alpha); g1_new(g_alpha);
    g1_mul(g_alpha, g1, alpha);
    pc_map(pp->gt, g_alpha, g2); 
}



void keygen_naive_oe(struct alp_pp_naive_oe pp, struct alp_sk_oe *sk, struct node *tree_root, bn_t alpha) { 
    init_secret_key_alp_oe(pp.bound, sk);
    std::vector<policy_coefficient> lsss_vector;
    lsss_vector = std::vector<policy_coefficient>(); 
    share_secret(tree_root, alpha, pp.order, lsss_vector, true); 
    for (auto it = lsss_vector.begin(); it != lsss_vector.end(); it++){
        struct alp_sk_attr_oe D;
        init_secret_key_attr_alp_oe(pp.bound, &D);
        size_t attr_index = it -> leaf_index-1;
        bn_t rho[pp.bound];
        bn_null(rho[0]); bn_new(rho[0]);
        bn_set_dig(rho[0], 1);
        for (size_t i = 1; i < pp.bound; i++) {
            bn_null(rho[i]); bn_new(rho[i]);
            bn_t i_read; bn_null(i_read); bn_new(i_read);
            bn_set_dig(rho[i], it -> leaf_index);
            bn_set_dig(i_read, i);
            bn_mxp_basic(rho[i], rho[i], i_read, pp.order); 
        }
        bn_t r_i; bn_null(r_i); bn_new(r_i);
        bn_rand_mod(r_i, pp.order);
        g2_t u_0_tmp; g2_null(u_0_tmp); g2_new(u_0_tmp);
        g2_mul(u_0_tmp, pp.U2[0], r_i);

        g2_mul(D.D1, pp.g2, it -> share);  
        g2_add(D.D1, D.D1, u_0_tmp);
        g2_mul(D.D2, pp.g2, r_i);
        for (int j = 0; j < pp.bound-1; j++){
            bn_t rho_i; bn_null(rho_i); bn_new(rho_i);
            bn_neg(rho_i, rho[j+1]);
            g2_null(D.K[j]); g2_new(D.K[j]);                
            g2_mul(D.K[j], pp.U2[1], rho_i);
            g2_add(D.K[j], D.K[j], pp.U2[j+2]);
            g2_mul(D.K[j], D.K[j], r_i);
        }
        sk->D[attr_index] = D; 
    }
}

void encrypt_naive_oe(struct alp_pp_naive_oe pp, bn_t *p_Coeffs, struct alp_ciphertext_oe *C) {
    gt_null(C->C0); gt_new(C->C0);
    g1_null(C->C1); g1_null(C->C1);
    g1_null(C->C2); g1_null(C->C2);
    g1_null(C->C3); g1_null(C->C3);

    bn_t s; bn_null(s); bn_new(s);
    bn_rand_mod(s, pp.order);
    //Naive imp, mul_sim here for optimisation
    gt_exp(C->C0, pp.gt, s);
    g1_mul(C->C1, pp.g1, s);
    g1_copy(C->C2, pp.U1[0]);
    g1_set_infty(C->C3);
    for (int i = 0; i < pp.bound; i++) {
        g1_t u_tmp; g1_null(u_tmp); g1_new(u_tmp);
        //cout << "pp.U1[" << i+1 << "]\n";
        //g1_print(pp.U1[i+1]);
        //cout << "p_Coeffs[" << i << "]\n";
        //bn_print(p_Coeffs[i]); 
        g1_mul(u_tmp, pp.U1[i+1], p_Coeffs[i]);
        g1_add(C->C2, C->C2, u_tmp);
        if (i < pp.bound-1) {
            g1_t h_tmp; g1_null(h_tmp); g1_new(h_tmp);
            g1_mul(h_tmp, pp.H1[i], p_Coeffs[i]);
            g1_add(C->C3, C->C3, h_tmp);
        }
        //cout << "u_tmp " << i << "\n";
        //g1_print(u_tmp);
    }

    g1_mul(C->C2, C->C2, s);
    g1_mul(C->C3, C->C3, s);


}


void coeff_array(bn_t *p_Coeffs, const bn_t *attributes, int bound){
    bn_t coeff_matrix[bound*bound];
    ind_coeff_bn(coeff_matrix, attributes, bound);
    for (size_t i = 0; i < bound; i++){
        bn_null(p_Coeffs[i]); bn_new(p_Coeffs[i])
        bn_copy(p_Coeffs[i], coeff_matrix[bound*(bound-1)+i]);
    }
    for (size_t i = 0; i < bound*bound; i++) {
        bn_free(coeff_matrix[i]);
    }
}

void decrypt_naive_oe(struct alp_pp_naive_oe pp, alp_sk_oe sk, alp_ciphertext_oe C, bn_t *attributes, struct node tree_root) {
    bn_t p_Coeffs[pp.bound];
    coeff_array(p_Coeffs, attributes, pp.bound);
    lsss_vector = std::vector<policy_coefficient>();
    lsss_vector = recover_coefficients(&tree_root, attributes, pp.bound-1);
    //gt_t share_points[pp.bound-1]; 
    gt_t result; gt_null(result); gt_new(result); gt_set_unity(result); 
    for (auto it = lsss_vector.begin(); it != lsss_vector.end(); it++) {
        size_t attr_index = it -> leaf_index-1;
        g2_t decrypt_d; g2_null(decrypt_d); g2_new(decrypt_d);
        g2_t Ky; g2_null(Ky); g2_new(Ky); g2_set_infty(Ky);
        for (int j = 0; j < pp.bound-1; j++){
            g2_t Ky_j; g2_null(Ky_j); g2_new(Ky_j);
            g2_mul(Ky_j, sk.D[attr_index].K[j], p_Coeffs[j+1]);
            g2_add(Ky, Ky, Ky_j);
        }
        g2_add(decrypt_d, sk.D[attr_index].D1, Ky);
        g2_t share_tmp2; 
        gt_t inv_tmp; gt_null(inv_tmp); gt_new(inv_tmp);
        gt_t share_point; gt_null(share_point); gt_new(share_point);

        pc_map(share_point, C.C1, decrypt_d);
        pc_map(inv_tmp, C.C2, sk.D[attr_index].D2);
        gt_inv(inv_tmp, inv_tmp);
        gt_mul(share_point, share_point, inv_tmp);
        gt_exp(share_point, share_point, it -> coeff);
        gt_mul(result, result, share_point);
    }  
    gt_inv(result, result);  
    gt_mul(result, result, C.C0);
    //cout << "res\n";
    int cmp  = gt_is_unity(result); 
    if (cmp != 1) {
        cout << "Value of result after decrypt\n";
        gt_print(result);
    } 
}

void decrypt_g_oe(struct alp_pp_naive_oe pp, alp_sk_oe sk, alp_ciphertext_oe C, bn_t *attributes, struct node tree_root) {
    bn_t p_Coeffs[pp.bound];
    coeff_array(p_Coeffs, attributes, pp.bound);
    lsss_vector = std::vector<policy_coefficient>();
    lsss_vector = recover_coefficients(&tree_root, attributes, pp.bound-1);
    //gt_t share_points[pp.bound-1]; 
    gt_t result; gt_null(result); gt_new(result); gt_set_unity(result); 
    for (auto it = lsss_vector.begin(); it != lsss_vector.end(); it++) {
        size_t attr_index = it -> leaf_index-1;
        g2_t decrypt_d; g2_null(decrypt_d); g2_new(decrypt_d);
        g2_t Ky; g2_null(Ky); g2_new(Ky); g2_set_infty(Ky);
        for (int j = 0; j < pp.bound-1; j++){
            g2_t Ky_j; g2_null(Ky_j); g2_new(Ky_j);
            g2_mul(Ky_j, sk.D[attr_index].K[j], p_Coeffs[j+1]);
            g2_add(Ky, Ky, Ky_j);
        }
        g2_add(decrypt_d, sk.D[attr_index].D1, Ky);
        g2_t share_tmp2; 
        gt_t inv_tmp; gt_null(inv_tmp); gt_new(inv_tmp);
        g1_t inv_C2; g1_null(inv_C2); g1_new(inv_C2);
        g1_t C1_share; g1_null(C1_share); g1_new(C1_share);

        gt_t share_point; gt_null(share_point); gt_new(share_point);
        g1_mul(C1_share, C.C1, it -> coeff);
       
        pc_map(share_point, C1_share, decrypt_d);
        g1_mul(inv_C2, C.C2, it -> coeff);
        g1_neg(inv_C2, inv_C2);
        pc_map(inv_tmp, inv_C2, sk.D[attr_index].D2);
        //gt_inv(inv_tmp, inv_tmp);
        gt_mul(share_point, share_point, inv_tmp);
        //gt_exp(share_point, share_point, it -> coeff);
        gt_mul(result, result, share_point);
    }  
    gt_inv(result, result);  
    gt_mul(result, result, C.C0);
    //cout << "res\n";
    int cmp  = gt_is_unity(result); 
    if (cmp != 1) {
        cout << "Value of result after decrypt\n";
        gt_print(result);
    } 
}

void print_secret_key_oe(struct alp_sk_oe sk, int bound) {
    for (int i = 0; i < bound-1; i++) {
        cout << "sk.D[" << i << "]" << ".D1\n"; 
        g2_print(sk.D[i].D1);
        cout << "sk.D[" << i << "]" << ".D2\n";
        g2_print(sk.D[i].D2);
        for (int j = 0; j < bound-1; j++){
            cout << "sk.D[" << i << "].K[" << j << "]\n";
            g2_print(sk.D[i].K[j]);
        } 
    }
}

void print_public_params(struct alp_pp_naive_oe pp, int bound) {
    cout << "pp.g1\n";
    g1_print(pp.g1);
    cout << "pp.g2\n";
    g2_print(pp.g2);
    cout << "pp.gt\n";
    gt_print(pp.gt);
    for (int i = 0; i < bound+1; i++) {
        cout << "pp.U1[" << i << "]\n";
        g1_print(pp.U1[i]);
        cout << "pp.U2[" << i << "]\n";
        g2_print(pp.U2[i]);
        if (i < bound) {
        cout << "pp.H1[" << i << "]\n";
            g1_print(pp.H1[i]);
            cout << "pp.H2[" << i << "]\n";
            g2_print(pp.H2[i]);
        }
    } 
} 

void print_share_component(struct alp_sk_oe sk, struct alp_pp_naive_oe pp, bn_t *p_coeffs, bn_t share, bn_t r, int bound, int attr_index, g2_t res) {
    gt_null(res); gt_new(res); g2_set_infty(res);
    //prod k_i_j ^ y_j
    //for (int i = 0; i < bound; i++) {
    //    g2_t u_y; gt_null(u_y); gt_new(u_y);
    //    g2_mul(u_y, pp.U2[i+1], p_coeffs[i]);
    //    g2_add(res, res, u_y)
    //}
    for (int i = 0; i < bound-1; i++) {
        g2_t k_y; gt_null(k_y); gt_new(k_y);
        g2_mul(k_y, sk.D[attr_index].K[i], p_coeffs[i+1]);
        g2_add(res, res, k_y);
    }
    //u0
    //g2_add(res, res, pp.U2[0]);
    //g2_mul(res, res, r);
    //u_0^r_i
    //g2_t g_share; gt_null(g_share); gt_new(g_share);
    //g2_mul(g_share, pp.g2, share);
    g2_add(res, res, sk.D[attr_index].D1);
    //g2_add(res, res, g_share);
    cout << "g^lamda*(prod u_j^y_j)^r_i\n";
    g2_print(res);
}