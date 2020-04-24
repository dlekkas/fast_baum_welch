#include <chrono>
#include <iostream>
#include <numeric>
#include <algorithm>

#include "../include/tsc_x86.h"
#include "../include/benchmark.h"


void perf_test_rdtscp(const std::string& init_file, std::vector<int>& observations,
		compute_func baum_welch, int n_runs, int n_iter, std::ostream& xout) {

	std::vector<double> cycles_list;

	init_tsc();

	for (auto i = 0; i < n_iter; i++) {
		HMM model(init_file);

		uint64_t start = start_tsc();
		for (auto j = 0; j < n_runs; j++) {
			baum_welch(model.M, model.N, observations.size(), observations.data(), \
				   model.pi.data(), model.A, model.B);
		}
		uint64_t end = stop_tsc();
		uint64_t cycles = (end - start) / (double) n_runs;

		cycles_list.emplace_back(cycles);
	}


	std::sort(cycles_list.begin(), cycles_list.end());

	double median = (n_iter % 2 == 0) ? cycles_list[n_iter/2] :
			(cycles_list[n_iter/2] + cycles_list[n_iter/2+1]) / 2;

	double average = std::accumulate(cycles_list.begin(),
			cycles_list.end(), 0.0) / cycles_list.size();

	xout << "[RDTSCP] Cycles elapsed (MIN): " << cycles_list.front() << std::endl;
	xout << "[RDTSCP] Cycles elapsed (MAX): " << cycles_list.back() << std::endl;
	xout << "[RDTSCP] Cycles elapsed (Mean): " << average << std::endl;
	xout << "[RDTSCP] Cycles elapsed (Median): " << median << std::endl;

}




void perf_test_chrono(const std::string& init_file, std::vector<int>& observations,
		compute_func baum_welch, int n_runs, int n_iter, std::ostream& xout) {

	std::vector<double> time_list;

	for (auto i = 0; i < n_iter; i++) {
		HMM model(init_file);
		auto begin = std::chrono::steady_clock::now();
		for (auto j = 0; j < n_runs; j++) {
			baum_welch(model.M, model.N, observations.size(), observations.data(), \
				   model.pi.data(), model.A, model.B);
		}
		auto end = std::chrono::steady_clock::now();

		auto duration_us = std::chrono::duration_cast
			<std::chrono::microseconds>(end - begin).count() / n_runs;

		time_list.emplace_back(duration_us);
	}

	std::sort(time_list.begin(), time_list.end());

	double median = (n_iter % 2 == 0) ? time_list[n_iter/2] :
			(time_list[n_iter/2] + time_list[n_iter/2+1]) / 2;

	double average = std::accumulate(time_list.begin(),
			time_list.end(), 0.0) / time_list.size();

	xout << "[CHRONO] Time (MIN): " << time_list.front() << " us" << std::endl;
	xout << "[CHRONO] Time (MAX): " << time_list.back() << " us" << std::endl;
	xout << "[CHRONO] Time (Mean): " << average << " us" << std::endl;
	xout << "[CHRONO] Time (Median): " << median << " us" << std::endl;

}
