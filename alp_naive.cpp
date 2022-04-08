#include <cstdio>
#include <string>
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
void setup_naive_oe(struct alp_pp_naive_oe *pp, bn_t order) {
    g1_t g1; g1_null(g1); g1_new(g1);
    g2_t g2; g2_null(g2); g2_new(g2);
    g2_null(g2); g2_new(g2);

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
        g1_mul(pp->U[i], g1, alpha_i);
        g1_mul(pp->H[i], g1, beta_i);
    }
    g1_copy(pp->g1, g1);
    g2_copy(pp->g2, g2);
    pc_map(pp->gt, g1, g2);
    bn_t alpha; bn_null(alpha); bn_new(alpha);
    bn_rand_mod(alpha, order);
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

    //Public params

    bn_t order;
    pc_param_set_any();
    pc_param_print();
    pc_get_ord(order);

    
    struct alp_pp_naive_oe pp;
    setup_naive_oe(&pp, order);

    gt_print(pp.gt);
}


int main (int argc, char **argv) {
    test(2);
    return 0;
}