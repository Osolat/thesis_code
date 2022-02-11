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

#include <stdint.h>


#include <relic/relic.h>

struct master_key {
    uint32_t N_ATTR;
    struct attribute * attributes;
    bn_t b;
    bn_t alpha;
};

struct attribute {
    uint32_t attr;
    bn_t b_attr;
};

int init_master_key(const uint32_t n_attr, struct master_key *m);

#endif


