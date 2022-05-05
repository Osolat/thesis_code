#include "t_function.h"

#include <iostream>

void lagrange_coeff(bn_t* res, size_t i, size_t N, bn_t x, bn_t order) {
    bn_set_dig(res[i - 1], 1);

    bn_t top_j;
    bn_null(top_j);
    bn_new(top_j);

    bn_t bot_i;
    bn_null(bot_i);
    bn_new(bot_i);
    bn_set_dig(bot_i, i);

    bn_t bot_j;
    bn_null(bot_j);
    bn_new(bot_j);

    for (size_t j = 0; j < N; j++) {
        if (j + 1 != i) {
            bn_set_dig(top_j, j + 1);
            bn_set_dig(bot_j, j + 1);
            bn_sub(top_j, x, top_j);
            bn_sub(bot_j, bot_i, bot_j);

            bn_mod_inv(bot_j, bot_j, order);

            bn_mul(bot_j, bot_j, top_j);

            bn_mod(res[i - 1], res[i - 1], order);

            bn_mul(res[i - 1], res[i - 1], bot_j);
        }
    }
    bn_mod(res[i - 1], res[i - 1], order);
}

void t_function_g2(g2_t* res, bn_t X, g2_t g2, g2_t* t_i_components, size_t n, bn_t order) {
    bn_t lags[n + 1];

    /*     bn_t* lags = (bn_t*)malloc((n + 1) * sizeof(bn_t));
     */
    for (size_t i = 0; i < n + 1; i++) {
        bn_null(lags[i]);
        bn_new(lags[i]);
        lagrange_coeff(lags, i + 1, n + 1, X, order);
    }
    g2_mul_sim_lot(*res, t_i_components, lags, n + 1);
    bn_t xn;
    bn_null(xn);
    bn_new(xn);
    g2_t temp;
    g2_null(temp);
    g2_new(temp);

    bn_set_dig(xn, n);
    bn_mxp(xn, X, xn, order);
    g2_mul(temp, g2, xn);
    g2_add(*res, temp, *res);
    /* for (size_t i = 0; i < n+1; i++)
    {
        bn_free(lags[i]);
    } */
}

void t_function_g2_pre(g2_t* res, bn_t X, g2_t* g2, g2_t* t_i_components, size_t n, bn_t order) {
    bn_t lags[n + 1];

    /*     bn_t* lags = (bn_t*)malloc((n + 1) * sizeof(bn_t));
     */
    for (size_t i = 0; i < n + 1; i++) {
        bn_null(lags[i]);
        bn_new(lags[i]);
        lagrange_coeff(lags, i + 1, n + 1, X, order);
    }
    g2_mul_sim_lot(*res, t_i_components, lags, n + 1);
    bn_t xn;
    bn_null(xn);
    bn_new(xn);
    g2_t temp;
    g2_null(temp);
    g2_new(temp);

    bn_set_dig(xn, n);
    bn_mxp(xn, X, xn, order);
    g2_mul_fix(temp, g2, xn);
    g2_add(*res, temp, *res);
    /* for (size_t i = 0; i < n+1; i++)
    {
        bn_free(lags[i]);
    } */
}