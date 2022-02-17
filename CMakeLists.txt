cmake_minimum_required(VERSION 3.21)
project(thesis_code)

set(CMAKE_CXX_STANDARD 17)

file(GLOB LEGACYFILES "/legacy/*.h" "/legacy/*.cpp")
file(GLOB LIBFILES "/lib/*.h" "/lib/*.cpp")

add_executable(thesis_code main.cpp ${LEGACYFILES} ${LIBFILES} bench_defs.h)
target_link_libraries(thesis_code /usr/local/lib/librelic_s.a)
target_link_libraries(thesis_code /usr/local/lib/libgmp.a)
