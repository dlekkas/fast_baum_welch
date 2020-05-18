#include <vector>
#include <tuple>
#include <iostream>

#include "../include/bw.h"
#include "../include/infra.h"

#define N_ITER 20
#define N_RUNS 1

using namespace std;

using InputParams= tuple<int, int, int>;
using Implementation = tuple<string, BaumWelch*>;

void run_benchmarks(const vector<InputParams> &inputs_var_M, const vector<Implementation> &implementations,
		const string &results_file_cycles, const string &results_file_time) {

	for (const auto& [M, N, S]: inputs_var_M) {
		std::cout << "M: " << M << ", N: " << N << ", S: " << S << std::endl;
		// Measurements regarding all the C-like implementations
		for (const auto [impl_tag, bw_func]: implementations) {
			perf_test_rdtscp(impl_tag, bw_func, M, N, S, N_RUNS, N_ITER, std::cout, true, results_file_cycles);
			perf_test_chrono(impl_tag, bw_func, M, N, S, N_RUNS, N_ITER, std::cout, true, results_file_time);
		}
	}
}

int main() {

	// entry: { <implementation-tag>, <baum-welch-function> }
	vector<Implementation> implementations {
		{"C++ Baseline", new BaumWelchCppBaseline()},
		{"C-like Baseline", new BaumWelchCBasic()},
		{"C Basic Opts", new BaumWelchCBasicOpts()},
		{"C More Opts", new BaumWelchCOptsV2()},
		{"C Loop Unroll", new BaumWelchCLoopUnroll()},
		{"C Basic Vectorized", new BaumWelchCVectBasic()}
	};

	// entry: { <n-states>, <n-emissions>, <observation-length> }

	vector<InputParams> inputs_var_M {
		{16,  64, 256},
		{32,  64, 256},
		{64,  64, 256},
		{128, 64, 256},
		{256, 64, 256},
		{512, 64, 256}
	};

	vector<InputParams> inputs_var_N {
		{128,  128, 256},
		{128,  256, 256},
		{128,  512, 256},
		{128, 1024, 256},
		{128, 2048, 256},
		{128, 4096, 256},
		{128, 8192, 256}
	};

	vector<InputParams> inputs_var_T {
		{128, 64,   64},
		{128, 64,  128},
		{128, 64,  256},
		{128, 64,  512},
		{128, 64, 1024},
		{128, 64, 2048},
		{128, 64, 4096}
	};

	run_benchmarks(inputs_var_M, implementations, "results_cycles_var_M.txt", "results_time_var_M.txt");
	run_benchmarks(inputs_var_N, implementations, "results_cycles_var_N.txt", "results_time_var_N.txt");
	run_benchmarks(inputs_var_T, implementations, "results_cycles_var_T.txt", "results_time_var_T.txt");

}
