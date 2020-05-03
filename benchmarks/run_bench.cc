#include <vector>
#include <tuple>
#include <iostream>

#include "../include/bw.h"
#include "../include/infra.h"

#define N_ITER 10
#define N_RUNS 1

using namespace std;

using InputParams= tuple<int,int,int>;
using Implementation = tuple<string, compute_func>;



int main() {

	// entry: { <implementation-tag>, <baum-welch-function> }
	vector<Implementation> implementations {
		{"C-like Baseline", &run_bw},
		{"Basic Opts", &run_bw_basic_opts}
	};

	// entry: : {<n-emissions>, <n-states>, <observation-length>}
	vector<InputParams> inputs {
		{16, 16, 256}, {16, 32, 256}, {64, 128, 256},
		{32, 32, 256}, {64, 64, 256}, {128, 128, 256},
		{32, 64, 256}, {64, 128, 256}, {128, 256, 256},
	};


	for (const auto& [N, M, S]: inputs) {
		// Measurements regarding the C++ baseline implementation
		perf_test_rdtscp("C++ Baseline", &baum_welch, N, M, S, N_RUNS, N_ITER, std::cout, true);
		perf_test_chrono("C++ Baseline", &baum_welch, N, M, S, N_RUNS, N_ITER, std::cout, true);

		// Measurements regarding all the C-like implementations
		for (const auto& [impl_tag, bw_func]: implementations) {
			perf_test_rdtscp(impl_tag, bw_func, N, M, S, N_RUNS, N_ITER, std::cout, true);
			perf_test_chrono(impl_tag, bw_func, N, M, S, N_RUNS, N_ITER, std::cout, true);
		}
	}

}
