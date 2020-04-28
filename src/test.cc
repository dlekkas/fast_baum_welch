#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>

#include "../include/bw.h"
#include "../include/benchmark.h"

#define N_ITERATIONS 1
#define N_RUNS 1


int main(int argc, char** argv) {
	std::cout << "Usage: ./baum_welch [<initialization file>] [<observations file>]" << std::endl;

	if (argc >= 3) {
		std::string init_file = argv[1];
		std::vector<int> observations;
		std::ifstream ifs(argv[2]);
		std::copy(std::istream_iterator<double>(ifs), std::istream_iterator<double>(),
				  std::back_inserter(observations));
		perf_test_chrono(init_file, observations, &run_bw, N_RUNS, N_ITERATIONS, std::cout);
	} else {
		perf_test_rdtscp_random(&baum_welch, 6, 2, N_RUNS, N_ITERATIONS, std::cout);
	}


	/* print all benchmarking results to a file */
	/*
	std::ofstream ofs("results.txt");
	perf_test_chrono(init_file, observations, &run_bw, N_RUNS, N_ITERATIONS, ofs);
	perf_test_rdtscp(init_file, observations, &run_bw, N_RUNS, N_ITERATIONS, ofs);
	*/

	/* print all benchmarking results to stdout */
	//perf_test_chrono(init_file, observations, &run_bw, N_RUNS, N_ITERATIONS, std::cout);
	//perf_test_rdtscp(init_file, observations, &run_bw, N_RUNS, N_ITERATIONS, std::cout);
}
