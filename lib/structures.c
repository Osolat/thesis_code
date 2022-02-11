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

#include <stdlib.h>
#include "structures.h"
#include "zp_arith.h"
#include "g1_arith.h"
#include "g2_arith.h"
#include "gt_arith.h"

int init_master_key(const uint32_t n_attr, struct master_key *m) {

    m->N_ATTR = n_attr;
    m->attributes = (struct attribute *) malloc(n_attr * sizeof(struct attribute));

    if (m->attributes == NULL) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

int init_public_key_cp(const uint32_t n_attr, struct public_key_cp *p) {

    p->N_ATTR = n_attr;
    p->attributes = (struct attribute_blind *) malloc (n_attr * sizeof(struct attribute_blind));

    if (p->attributes == NULL) {
        return EXIT_FAILURE;
    } else {

        g1_ar_init(&p->g);
        g1_ar_init(&p->B);
        g2_ar_init(&p->h);

        for (uint32_t i = 0; i < n_attr; i++)
            g1_ar_init(&p->attributes[i].g_b_attr);

        return EXIT_SUCCESS;
    }
}








