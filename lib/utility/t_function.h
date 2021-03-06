#include <stdint.h>

extern "C" {
#include <relic/relic.h>
}

void t_function_g2(g2_t* res, bn_t X, g2_t g2, g2_t* t_i_components, size_t n, bn_t order);
void t_function_g2_gap(g2_t* res, bn_t X, g2_t* g2, g2_t* t_i_components, size_t n, bn_t order);
void t_function_g1_gap(g1_t* res, bn_t X, g1_t* g2, g1_t* t_i_components, size_t n, bn_t order);