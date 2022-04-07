#include <stdint.h>

extern "C" {
#include <relic/relic.h>
}

void t_function_g2(g2_t* res, bn_t X, g2_t g2, g2_t* t_i_components, size_t n, bn_t order);