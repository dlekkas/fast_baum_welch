cmake_minimum_required(VERSION 3.10)
project(baum_welch_asl)

# Compilation flags on steroids
#set(CMAKE_CXX_FLAGS "-Wall -std=c++17 -mavx2 -mfma -O3 -flto -funroll-loops -march=native -mavx2 -mfma")

# Basic compilation flags - development
set(CMAKE_CXX_FLAGS "-Wall -std=c++17 -mavx2 -mfma -O3")

#set(CMAKE_CXX_FLAGS "-Wall -std=c++17 -Ofast -flto -funroll-loops -frename-registers -march=native -mavx2 -mfma")
#set(CMAKE_CXX_FLAGS "-Wall -O3 -std=c++17 -mavx2 -mfma -L /home/andreas/ETH/Advanced_Systems_Lab/pcm/pcm.so -lpcm")


# GNU COMPILER - GCC
#set(CMAKE_CXX_FLAGS "-Wall -std=c++17 -fprofile-generate -O3 -mavx2 -mfma")
#set(CMAKE_CXX_FLAGS "-Wall -std=c++17 -fprofile-use -O3 -mavx2 -mfma")

# INTEL COMPILER - ICC
#set(CMAKE_CXX_FLAGS "-Wall -std=c++17 -prof-gen -prof-dir=prof/ -O3 -mavx2 -mfma")
#set(CMAKE_CXX_FLAGS "-Wall -std=c++17 -prof-use -prof-dir=prof/ -O3 -mavx2 -mfma")



set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(INCLUDE_DIR "include")
include_directories(${INCLUDE_DIR})


file(GLOB IMPLEMENTATIONS "src/baum_welch/*.cc")
file(GLOB BENCHMARKS "benchmarks/*.cc")
file(GLOB SOURCES "src/*.cc")


foreach( bench_file ${BENCHMARKS} )

	string(REPLACE ".cc" "" exec_name ${bench_file})
	get_filename_component(exec_name ${exec_name} NAME)
	add_executable(${exec_name} ${SOURCES} ${IMPLEMENTATIONS} ${bench_file})

endforeach( bench_file ${BENCHMARKS} )
