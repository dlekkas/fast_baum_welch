#include <tuple>
#include <vector>
#include <iostream>

#include "../include/bw.h"
#include "../include/infra.h"
#include "../include/hmm.h"
#include "../include/generator.h"

#define N_ITERATIONS 10
#define N_RUNS 1

#define SEQ_LEN 128
#define M 64
#define N 64


using namespace std;
using Implementation = tuple<string, compute_func>;


int main() {

	// entry: { <implementation-tag> <baum-welch-function> }
	vector<Implementation> implementations {
		{"C-like Baseline", &run_bw},
		{"Basic Opts", &run_bw_basic_opts}
	};


	for (const auto& [impl_tag, bw_func]: implementations) {
		if (!IsValidImpl(bw_func)) {
			cout << "[" << impl_tag << "] Invalid implementation!" << endl;
		}
		perf_test_rdtscp(impl_tag, bw_func, M, N, SEQ_LEN, N_RUNS, N_ITERATIONS, std::cout);
		perf_test_chrono(impl_tag, bw_func, M, N, SEQ_LEN, N_RUNS, N_ITERATIONS, std::cout);
	}

}
