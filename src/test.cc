#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>

#include "../include/bw.h"
#include "../include/benchmark.h"

#define N_ITERATIONS 10
#define N_RUNS 1


int main(int argc, char** argv) {
    if (argc < 3) {
		std::cout << "Usage: ./baum_welch <initialization file> <observations file>" << std::endl;
		std::exit(1);
    }
	std::string init_file {argv[1]};

    std::vector<int> observations;
    std::ifstream ifs(argv[2]);
	std::copy(std::istream_iterator<double>(ifs), std::istream_iterator<double>(),
			  std::back_inserter(observations));


	/* print all benchmarking results to a file */
	/*
	std::ofstream ofs("results.txt");
	perf_test_chrono(init_file, observations, &run_bw, N_RUNS, N_ITERATIONS, ofs);
	perf_test_rdtscp(init_file, observations, &run_bw, N_RUNS, N_ITERATIONS, ofs);
	*/

	/* print all benchmarking results to stdout */
	perf_test_chrono(init_file, observations, &run_bw, N_RUNS, N_ITERATIONS, std::cout);
	perf_test_rdtscp(init_file, observations, &run_bw, N_RUNS, N_ITERATIONS, std::cout);
}
