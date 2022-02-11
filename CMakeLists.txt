cmake_minimum_required(VERSION 3.21)
project(thesis_code)

set(CMAKE_CXX_STANDARD 17)
add_executable(thesis_code main.cpp lib/structures.c lib/structures.h)
target_link_libraries(thesis_code /usr/local/lib/librelic_s_type3.a)
target_link_libraries(thesis_code /usr/local/lib/libgmp.a)
