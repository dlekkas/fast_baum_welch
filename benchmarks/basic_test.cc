#include <tuple>
#include <vector>
#include <iostream>

#include "../include/infra.h"
#include "../include/baum_welch.h"

#define N_ITERATIONS 25
#define N_RUNS 1

#define SEQ_LEN 128
#define M 64
#define N 64


using namespace std;
using Implementation = tuple<string, BaumWelch*>;


int main() {

	vector<Implementation> implementations {
		//{"C++ baseline", new BaumWelchCppBaseline()},
		//{"C++ opts", new BaumWelchCppOpts()},
<<<<<<< HEAD
		{"C basic", new BaumWelchCBasic()},
		{"C basic opts", new BaumWelchCBasicOpts()},
		{"C more opts", new BaumWelchCOptsV2()},
		{"C loop unrolling v1", new BaumWelchCLoopUnroll0()},
=======
		//{"C basic", new BaumWelchCBasic()},
		//{"C basic opts", new BaumWelchCBasicOpts()},
		//{"C more opts", new BaumWelchCOptsV2()},
		//{"C loop unrolling v1", new BaumWelchCLoopUnroll0()},
>>>>>>> f0ecfcf8779339dfacfd85a4b6336fe7ca60287e
		{"C loop unrolling opt", new BaumWelchCLoopUnroll()},
		{"C vectorized opt", new BaumWelchCVectOpt()},
		{"C vectorized v2", new BaumWelchCVectDim()},
	};

	for (auto [impl_tag, bw_func]: implementations) {
		if (!IsValidImpl(bw_func)) {
			cout << "[" << impl_tag << "] Invalid implementation!" << endl;
		}
		perf_test_rdtscp(impl_tag, bw_func, M, N, SEQ_LEN, N_RUNS, N_ITERATIONS, std::cout);
		perf_test_chrono(impl_tag, bw_func, M, N, SEQ_LEN, N_RUNS, N_ITERATIONS, std::cout);
	}

}
