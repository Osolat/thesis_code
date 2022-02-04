cmake_minimum_required(VERSION 3.21)
project(thesis_code C)

set(CMAKE_C_STANDARD 17)
add_executable(thesis_code main.c)
target_link_libraries(thesis_code /usr/local/lib/librelic_s_type3.a)
target_link_libraries(thesis_code /usr/local/lib/libgmp.a)
