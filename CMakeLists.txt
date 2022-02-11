cmake_minimum_required(VERSION 3.21)
project(thesis_code)

set(CMAKE_CXX_STANDARD 17)
add_executable(thesis_code main.cpp lib/structures.c lib/structures.h lib/gt_arith.c lib/gt_arith.h lib/g1_arith.c lib/g1_arith.h lib/g2_arith.h lib/g2_arith.c lib/pairing_arith.h lib/pairing_arith.c lib/zp_arith.c lib/zp_arith.h)
target_link_libraries(thesis_code /usr/local/lib/librelic_s_type3.a)
target_link_libraries(thesis_code /usr/local/lib/libgmp.a)
