#include <chrono>
#include <iostream>

#include "../include/tsc_x86.h"
#include "../include/benchmark.h"


void perf_test_rdtscp(HMM& model, std::vector<int>& observations, compute_func baum_welch) {
	init_tsc();

	uint64_t start = start_tsc();
    baum_welch(model.M, model.N, observations.size(), observations.data(), \
		   model.pi.data(), model.A, model.B);
	uint64_t end = stop_tsc();

	uint64_t cycles = end - start;
	std::cout << "[RDTSCP] Total cycles elapsed = " << cycles << std::endl;

}


void perf_test_chrono(HMM& model, std::vector<int>& observations, compute_func baum_welch) {

	auto begin = std::chrono::steady_clock::now();
    baum_welch(model.M, model.N, observations.size(), observations.data(), \
		   model.pi.data(), model.A, model.B);
	auto end = std::chrono::steady_clock::now();

	auto duration_us = std::chrono::duration_cast
		<std::chrono::microseconds>(end - begin).count();

	auto duration_ns = std::chrono::duration_cast
		<std::chrono::nanoseconds> (end - begin).count();

	std::cout << "[CHRONO] Time : " << duration_us << " (us)" << std::endl;
	std::cout << "[CHRONO] Time : " << duration_ns << " (ns)" <<  std::endl;
}
