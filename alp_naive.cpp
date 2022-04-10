#include <cstdio>
#include <string>
#include <cmath>
#include "bench_defs.h"
using namespace std;

int N_ATTR;

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
    for (auto i = 0; i < n; i++){
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

unsigned long long t[NTESTS];

void setup_naive_oe(struct alp_pp_naive_oe *pp, bn_t alpha, bn_t order) {
    g1_t g1; g1_null(g1); g1_new(g1);
    g2_t g2; g2_null(g2); g2_new(g2);

    g1_get_gen(g1);
    //cout << "g1\n";
    //g1_print(g1);
    g2_get_gen(g2);
    //cout << "g2\n";
    //g2_print(g2);
    init_public_params_alp_oe(N_ATTR, pp);     
    for(int i = 0; i < N_ATTR; i++) {
        bn_t alpha_i; bn_null(alpha_i); bn_new(alpha_i);
        bn_t beta_i; bn_null(beta_i); bn_new(beta_i);
        bn_rand_mod(alpha_i, order); bn_rand_mod(beta_i, order);
        g1_mul(pp->U1[i], g1, alpha_i);
        g2_mul(pp->U2[i], g2, alpha_i);
        g1_mul(pp->H1[i], g1, beta_i);
        g2_mul(pp->H2[i], g2, beta_i);
    }
    g1_copy(pp->g1, g1);
    g2_copy(pp->g2, g2);
    pc_map(pp->gt, g1, g2);
    gt_exp(pp->gt, pp->gt, alpha);
}

void test(int N_ATTR) {
    N_ATTR = N_ATTR;
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

    cout << "Public params\n";

    bn_t order;
    pc_param_set_any();
    pc_param_print();
    pc_get_ord(order);

    
    struct alp_pp_naive_oe pp;
    bn_t alpha; bn_null(alpha); bn_new(alpha);
    bn_rand_mod(alpha, order);
    setup_naive_oe(&pp, alpha, order);
   
 
 
    cout << "Key Gen\n";
    struct alp_sk_oe sk;
    struct node tree_root;
    vector<policy_coefficient> vector;
    tree_from_string(and_tree_formula(N_ATTR), &tree_root);
    init_secret_key_alp_oe(N_ATTR, &sk);
    for (size_t j = 0; j < 1; j++) {
        vector = std::vector<policy_coefficient>(); 
        share_secret(&tree_root, alpha, order, vector, true); 
        for (auto it = vector.begin(); it != vector.end(); it++){
            struct alp_sk_attr_oe D;
            init_secret_key_attr_alp_oe(N_ATTR, &D);
            size_t attr_index = it -> leaf_index-1;
            //cout << "leaf_index: " << it -> leaf_index << "\n";
            bn_t rho[N_ATTR];
            for (size_t i = 0; i < N_ATTR; i++) {
                bn_null(rho[i]); bn_new(rho[i]);
                bn_t i_read; bn_null(i_read); bn_new(i_read);
                bn_set_dig(rho[i], it -> leaf_index);
                bn_set_dig(i_read, i);
                bn_mxp_basic(rho[i], rho[i], i_read, order); 
                //cout << "rho[" << i << "]\n";
                //bn_print(rho[i]);
            }
            bn_t r_i; bn_null(r_i); bn_new(r_i);
            bn_rand_mod(r_i, order);
            g2_mul(D.D1, pp.g2, it -> share); 
            g2_t u_0_tmp; g2_null(u_0_tmp); g2_new(u_0_tmp);
            g2_mul(u_0_tmp, pp.U2[0], r_i);
            g2_add(D.D1, D.D1, u_0_tmp);
            //cout << "D_" << it -> leaf_index << "_1\n";
            //g2_print(D.D1);

            g2_mul(D.D2, pp.g2, it -> share);
            //cout << "D_" << it -> leaf_index << "_2\n";
            //g2_print(D.D2);

            for (int j = 0; j < N_ATTR-1; j++){
                bn_t rho_i; bn_null(rho_i); bn_new(rho_i);
                //1/rho_{i,1}
                bn_mod_inv(rho_i, rho[0], order);
                //rho_{i,j+1}/rho_{i,1}
                bn_mul(rho_i, rho_i, rho[j+1]);
                //-rho_{i,j+1}/rho_{i,1}
                bn_neg(rho_i, rho_i);
                bn_mod(rho_i, rho_i, order);
                //cout << "rho[" << attr_index << "][" << j << "]\n";

                g2_null(D.K[(N_ATTR-1)*attr_index+j]); g2_new(D.K[(N_ATTR-1)*attr_index+j]);
                //bn_print(rho_i);
                //u_1^(-rho_{i,j+1}/rho_{i,1})
                g2_mul(D.K[(N_ATTR-1)*attr_index+j], pp.U2[0], rho_i);
                //cout << "sheesh\n";
                //u_1^(-rho_{i,j+1}/rho_{i,1})*u_2
                g2_add(D.K[(N_ATTR-1)*attr_index+j], D.K[(N_ATTR-1)*attr_index+j], pp.U2[1]);
                //(u_1^(-rho_{i,j+1}/rho_{i,1})*u_2)^r_i
                g2_mul(D.K[(N_ATTR-1)*attr_index+j], D.K[(N_ATTR-1)*attr_index+j],r_i);
                //cout << "K[" << attr_index << "][" << j << "]\n";
                //g2_print(D.K[(N_ATTR-1)*attr_index+j]);
            }
            sk.D[attr_index] = D; 
        }
    } 

    cout << "Encrypt\n";
    bn_t attributes[N_ATTR];
    std::vector<int> attr_vector(N_ATTR);
    bn_t p_Coeffs[N_ATTR];
    for (size_t i = 0; i < N_ATTR; i++) {
        bn_null(attributes[i]); 
        bn_new(attributes[i]); 
        bn_set_dig(attributes[i], i + 1);
        attr_vector[i] = i+1;
    }
    //TODO: This is the q+1 = n_attr case, if we don't encrypt under all attributes, this does not work.
    std::vector<std::vector<int>> coeff_vector(N_ATTR+1, std::vector<int> (N_ATTR+1, 0));
    cout << "yeesh\n";
    ind_coeff(coeff_vector, attr_vector, N_ATTR);
    for (size_t i = 0; i < N_ATTR; i++){
        bn_null(p_Coeffs[i]; bn_new(p_Coeffs[i]));
        int y = coeff_vector[N_ATTR-1][i];
        //cout << "y_" << i << ": " << y <<  "\n";
        bn_set_dig(p_Coeffs[i], y);
    }
    bn_t s; bn_null(s); bn_new(s);
    gt_t C0; gt_null(C0); gt_new(C0);
    g1_t C1; g1_null(C1); g1_new(C1);
    g1_t C2; g2_null(C2); g2_null(C2);
    g1_t C3; g2_null(C3); g2_null(C3);


    bn_rand_mod(s, order);
    //Naive imp, mul_sim here for optimisation
    gt_exp(C0, pp.gt, s);
    g1_mul(C1, pp.g1, s);
    g1_copy(C2, pp.U1[0]);
    g1_copy(C3, pp.H1[0]);
    for (int i = 1; i < N_ATTR; i++) {
        g1_t u_tmp; g1_null(u_tmp); g1_new(u_tmp);
        g1_t h_tmp; g1_null(h_tmp); g1_new(h_tmp);
        g1_mul(u_tmp, pp.U1[i], p_Coeffs[i]);
        g1_mul(h_tmp, pp.H1[i], p_Coeffs[i]);
        g1_add(C2, C2, u_tmp);
        g1_add(C3, C3, h_tmp);
    }
    g1_mul(C2, C2, s);
    g1_mul(C3, C3, s);

    cout << "Decrypt\n";
    try {
        check_satisfiability(&tree_root, attributes, N_ATTR);
        //std::cout << "Satisfiable with correct attributes" << std::endl;
    } catch (struct TreeUnsatisfiableException *e) {
        std::cout << e->what() << std::endl;
    }

    vector = std::vector<policy_coefficient>();
    vector = recover_coefficients(&tree_root, attributes, N_ATTR);
    for (auto it = vector.begin(); it != vector.end(); it++) {
        size_t attr_index = it -> leaf_index-1;
        g2_t decrypt_d; g2_null(decrypt_d); g2_new(decrypt_d);
        g2_copy(decrypt_d, sk.D->D1);
        for (int j = 1; j < N_ATTR; j++){
            g2_t Ky; g2_null(Ky); g2_new(Ky);
            //g2_mul(Ky, sk.D[attr_index].K[attr_index][j], p_Coeffs[j]);
        }
    }    
}


int main (int argc, char **argv) {
    test(8);
    return 0;
}