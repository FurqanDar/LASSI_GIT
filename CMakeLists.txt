cmake_minimum_required(VERSION 3.10)
project(LASSI_SRC C)

set(CMAKE_CXX_COMPILER "clang")
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

include_directories(.)

add_executable(lassi
        src/cluster.c
        src/cluster.h
        src/energy.c
        src/energy.h
        src/global.h
        src/initialize.c
        src/initialize.h
        src/main.c
        src/mcmove.c
        src/mcmove.h
        src/parsekey.c
        src/parsekey.h
        src/print.c
        src/print.h
        src/structure.c
        src/structure.h)

target_link_libraries(lassi m)