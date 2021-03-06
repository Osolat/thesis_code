
#ifndef BENCH_DEFS_H
#define BENCH_DEFS_H

#include <assert.h>
#include <sched.h>
#include <stdbool.h>
#include <unistd.h>

#include "lib/policy/policy_tree.h"
#include "lib/zp_arith.h"
#include "lib/structures.h"
#include "lib/pairing_arith.h"
#include "lib/gt_arith.h"
#include "lib/g1_arith.h"
#include "lib/utility/t_function.h"

#include "legacy/l_zattributelist.h"
#include "legacy/l_zconstants.h"
#include "legacy/l_zinteger.h"
#include "legacy/l_zpolicy.h"
#include "legacy/l_zbytestring.h"
#include "legacy/l_zfunctioninput.h"
#include "legacy/l_zobject.h"

#include "legacy/arith/l_zgroup.h"
#include "legacy/arith/l_zelement_bp.h"
#include "legacy/lsss/l_zlsss.h"

extern "C" {
#include "legacy/arith/l_zelement.h"
}

/* Number of tests per operation during measurement */

#define NTESTS 2000
#define g1_pre_sim_switchpoint 32
#define g2_pre_sim_switchpoint 32
#define exponent_bits_exceed_breakpoint 47
#endif
