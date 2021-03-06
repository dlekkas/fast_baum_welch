#include <tuple>
#include <vector>
#include <iostream>

#include "../include/infra.h"
#include "../include/baum_welch.h"

#define N_ITERATIONS 20
#define N_RUNS 1

#define SEQ_LEN 256
#define M 128
#define N 128


using namespace std;
using Implementation = tuple<string, BaumWelch*>;


int main() {

	vector<Implementation> implementations {
		//{"C++ baseline", new BaumWelchCppBaseline()},
		//{"C++ opts", new BaumWelchCppOpts()},
		//{"C basic", new BaumWelchCBasic()},
		//{"C basic opts", new BaumWelchCBasicOpts()},
		//{"C more opts", new BaumWelchCOptsV2()},
		//{"C Manos", new BaumWelchCOptsManos()},
		//{"C Blocking", new BaumWelchCOptsBlocking()},
		{"C loop unrolling opt", new BaumWelchCLoopUnroll()},
		{"C loop unrolling + Blocking", new BaumWelchCBlocking_Unroll()},
		//{"C vectorized opt", new BaumWelchCVectOpt()},
		{"C vectorized v3", new BaumWelchCVectDim2()},
		{"C vectorized+blocking", new BaumWelchCVectBlocking()},
		//{"C vectorized with unroll", new BaumWelchCVectUnroll()},
	};

	for (auto [impl_tag, bw_func]: implementations) {
		if (!IsValidImpl(bw_func)) {
			cout << "[" << impl_tag << "] Invalid implementation!" << endl;
		}
		perf_test_rdtscp(impl_tag, bw_func, M, N, SEQ_LEN, N_RUNS, N_ITERATIONS, std::cout);
		perf_test_chrono(impl_tag, bw_func, M, N, SEQ_LEN, N_RUNS, N_ITERATIONS, std::cout);
	}

}
