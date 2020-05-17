cmake_minimum_required(VERSION 3.10)
project(baum_welch_asl)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-Wall -O3 -mavx2 -mfma")


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
