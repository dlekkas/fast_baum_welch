#include <tuple>
#include <vector>
#include <iostream>

#include "../include/bw.h"
#include "../include/infra.h"
#include "../include/hmm.h"
#include "../include/generator.h"

#define N_ITERATIONS 25
#define N_RUNS 1

#define SEQ_LEN 128
#define M 64
#define N 64


using namespace std;
using Implementation = tuple<string, compute_func>;
using Implementation_c = tuple<string, compute_func1>;
using Implementation_cpp = tuple<string, compute_func2>;



int main() {

	vector<Implementation_cpp> implementations_new {
		{"C++ Baseline", &baum_welch},
		{"C++ Basic Opts", &baum_welch_opts}
	};

	for (const auto& [impl_tag, bw_func]: implementations_new) {
		if (!IsValidImpl(bw_func)) {
			cout << "[" << impl_tag << "] Invalid implementation!" << endl;
		}
		perf_test_rdtscp(impl_tag, bw_func, M, N, SEQ_LEN, N_RUNS, N_ITERATIONS, std::cout);
		perf_test_chrono(impl_tag, bw_func, M, N, SEQ_LEN, N_RUNS, N_ITERATIONS, std::cout);
	}

	vector<Implementation> implementations {
		{"C-like Baseline", &run_bw},
		{"C-like Basic Opts", &run_bw_basic_opts},
		{"C-like More Opts", &run_bw_opts_v2},
	};


	for (const auto& [impl_tag, bw_func]: implementations) {
		if (!IsValidImpl(bw_func)) {
			cout << "[" << impl_tag << "] Invalid implementation!" << endl;
		}
		perf_test_rdtscp(impl_tag, bw_func, M, N, SEQ_LEN, N_RUNS, N_ITERATIONS, std::cout);
		perf_test_chrono(impl_tag, bw_func, M, N, SEQ_LEN, N_RUNS, N_ITERATIONS, std::cout);
	}


	vector<Implementation_c> implementations_opt {
		{"C-like Loop Unroll", &bw_loop_unroll},
		{"C-like Loop Unroll Opts", &bw_loop_unroll_opt},

	};


	for (const auto& [impl_tag, bw_func]: implementations_opt) {
		if (!IsValidImpl(bw_func)) {
			cout << "[" << impl_tag << "] Invalid implementation!" << endl;
		}
		perf_test_rdtscp(impl_tag, bw_func, M, N, SEQ_LEN, N_RUNS, N_ITERATIONS, std::cout);
		perf_test_chrono(impl_tag, bw_func, M, N, SEQ_LEN, N_RUNS, N_ITERATIONS, std::cout);
	}

}
