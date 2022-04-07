#include <stdint.h>

extern "C" {
#include <relic/relic.h>
}

/**
 * Evaluates lagrange basis polynomial function \Delta_i_N (X)
 * @param[in] res		- pointer to res value that will have final result
 * @param[in] i			- pointer to res value that will have final result
 * @param[in] N			- pointer to res value that will have final result
 * @param[in] x			- pointer to res value that will have final result
 */
void lagrange_coeff(bn_t* res, size_t i, size_t N, bn_t x, bn_t order);

void t_function_g2(g2_t* res, bn_t X, g2_t g2, g2_t* t_i_components, size_t n, bn_t order);