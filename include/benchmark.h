#ifndef BENCHMARK_H
#define BENCHMARK_H

#include "hmm.h"

using compute_func = void (*)(int, int, int, int*, double*, double**, double**);

void perf_test_rdtscp(HMM& model, std::vector<int>& observations, compute_func baum_welch);

void perf_test_chrono(HMM& model, std::vector<int>& observations, compute_func baum_welch);

#endif
