#include <vector>
#include <tuple>
#include <iostream>

#include "../include/bw.h"
#include "../include/infra.h"

#define N_ITER 10
#define N_RUNS 1

using namespace std;

using InputParams= tuple<int, int, int>;
using Implementation = tuple<string, compute_func>;



int main() {

	// entry: { <implementation-tag>, <baum-welch-function> }
	vector<Implementation> implementations {
		{"C-like Baseline", &run_bw},
		{"Basic Opts", &run_bw_basic_opts}
	};

	// entry: : {<n-states>, <n-emissions>, <observation-length>}
	/*
	vector<InputParams> inputs {
		{16, 16, 256}, {16, 32, 256},
		{32, 32, 256}, {32, 64, 256},
		{64, 64, 256}, {64, 128, 256},
		{128, 128, 256}, {128, 256, 256}
	};
	*/
	vector<InputParams> inputs {
		{16, 16, 256}, {32, 32, 256}, {64, 64, 256},
		{128, 128, 256}, {256, 256, 256}, {512, 512, 256}
	};

	string results_file_cycles = "results_cycles.txt";
	string results_file_time   = "results_time.txt";
	for (const auto& [M, N, S]: inputs) {
		std::cout << "M: " << M << ", N: " << N << ", S: " << S << std::endl;
		// Measurements regarding the C++ baseline implementation
		perf_test_rdtscp("C++ Baseline", &baum_welch, M, N, S, N_RUNS, N_ITER, std::cout, true, results_file_cycles);
		perf_test_chrono("C++ Baseline", &baum_welch, M, N, S, N_RUNS, N_ITER, std::cout, true, results_file_time);

		// Measurements regarding all the C-like implementations
		for (const auto& [impl_tag, bw_func]: implementations) {
			perf_test_rdtscp(impl_tag, bw_func, M, N, S, N_RUNS, N_ITER, std::cout, true, results_file_cycles);
			perf_test_chrono(impl_tag, bw_func, M, N, S, N_RUNS, N_ITER, std::cout, true, results_file_time);
		}

		// Measurements regarding the C++ baseline implementation
		perf_test_rdtscp("C++ Baseline Opts", &baum_welch_opts, M, N, S, N_RUNS, N_ITER, std::cout, true,
			results_file_cycles);
		perf_test_chrono("C++ Baseline Opts", &baum_welch_opts, M, N, S, N_RUNS, N_ITER, std::cout, true,
			results_file_time);
	}

}
