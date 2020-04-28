#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>

#include "../include/bw.h"
#include "../include/benchmark.h"

#define N_ITERATIONS 10
#define N_RUNS 1


int main(int argc, char** argv) {

	/*
	 *  This is how we would get initialization parameters + observations from file
		std::string init_file { "init_params.txt" };
		std::string obs_file { "observations.txt" };
		perf_test_chrono(init_file, obs_file, &run_bw, N_RUNS, N_ITERATIONS, std::cout);
		perf_test_rdtscp(init_file, obs_file, &run_bw, N_RUNS, N_ITERATIONS, std::cout);
	*/

	int n_states = 5;
	int n_emissions = 3;

	// Measurements regarding the C++ baseline implementation
	perf_test_rdtscp("Baseline implementation", &baum_welch,
			n_states, n_emissions, N_RUNS, N_ITERATIONS, std::cout);
	perf_test_chrono("Baseline implementation", &baum_welch,
			n_states, n_emissions, N_RUNS, N_ITERATIONS, std::cout);

	// Measurements regarding the C-like implementation
	perf_test_rdtscp("Basic opts", &run_bw, n_states, n_emissions, N_RUNS, N_ITERATIONS, std::cout);
	perf_test_chrono("Basic opts", &run_bw, n_states, n_emissions, N_RUNS, N_ITERATIONS, std::cout);

}
