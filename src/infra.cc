#include <chrono>
#include <iostream>
#include <fstream>
#include <numeric>
#include <iterator>
#include <algorithm>

#include "../include/tsc_x86.h"
#include "../include/benchmark.h"
#include "../include/infra.h"

#define SEQ_LEN 50


void perf_test_rdtscp(const std::string& impl_tag, const std::string& init_file,
		const std::string& obs_file, compute_func baum_welch, int n_runs,
		int n_iter, std::ostream& xout) {

	std::vector<double> cycles_list;
	std::vector<int> observations;
	std::ifstream ifs(obs_file);
	std::copy(std::istream_iterator<int>(ifs), std::istream_iterator<int>(),
			  std::back_inserter(observations));

	init_tsc();

	HMM base_model(init_file);
	for (auto i = 0; i < n_iter; i++) {
		HMM model(base_model);

		uint64_t start = start_tsc();
		for (auto j = 0; j < n_runs; j++) {
			baum_welch(model.M, model.N, observations.size(), observations.data(), \
				   model.pi.data(), model.A, model.B);
		}
		uint64_t end = stop_tsc();
		uint64_t cycles = (end - start) / (double) n_runs;

		cycles_list.emplace_back(cycles);
	}

	Benchmark bench(cycles_list, impl_tag, "cycles", base_model.N, base_model.M, observations.size());
	bench.BeautyPrint(xout);

}


void perf_test_rdtscp(const std::string& impl_tag, compute_func baum_welch,
		int M, int N, int n_runs, int n_iter, std::ostream& xout) {

	std::vector<int> observations = uniform_emission_sample(SEQ_LEN, N);
	std::vector<double> cycles_list;

	init_tsc();

	HMM base_model(M, N);
	for (auto i = 0; i < n_iter; i++) {

		HMM model(base_model);

		uint64_t start = start_tsc();
		for (auto j = 0; j < n_runs; j++) {
			baum_welch(model.M, model.N, observations.size(), observations.data(), \
				   model.pi.data(), model.A, model.B);
		}
		uint64_t end = stop_tsc();
		uint64_t cycles = (end - start) / (double) n_runs;

		cycles_list.emplace_back(cycles);
	}


	Benchmark bench(cycles_list, impl_tag, "cycles", N, M, observations.size());
	bench.BeautyPrint(xout);
}


void perf_test_rdtscp(const std::string& impl_tag, compute_func2 baum_welch,
		int M, int N, int n_runs, int n_iter, std::ostream& xout) {

	std::vector<int> observations = uniform_emission_sample(SEQ_LEN, N);
	std::vector<double> cycles_list;

	init_tsc();

	HMM base_model(M, N);
	for (auto i = 0; i < n_iter; i++) {

		HMM model(base_model);

		uint64_t start = start_tsc();
		for (auto j = 0; j < n_runs; j++) {
			baum_welch(model.transition, model.emission, model.pi, observations);
		}
		uint64_t end = stop_tsc();
		uint64_t cycles = (end - start) / (double) n_runs;

		cycles_list.emplace_back(cycles);
	}

	Benchmark bench(cycles_list, impl_tag, "cycles", N, M, observations.size());
	bench.BeautyPrint(xout);

}


void perf_test_chrono(const std::string& impl_tag, const std::string& init_file,
		const std::string& obs_file, compute_func baum_welch, int n_runs,
		int n_iter, std::ostream& xout) {

	std::vector<double> time_list;

	std::vector<int> observations;
	std::ifstream ifs(obs_file);
	std::copy(std::istream_iterator<int>(ifs), std::istream_iterator<int>(),
			  std::back_inserter(observations));

	HMM base_model(init_file);
	for (auto i = 0; i < n_iter; i++) {
		HMM model(base_model);

		auto begin = std::chrono::steady_clock::now();
		for (auto j = 0; j < n_runs; j++) {
			baum_welch(model.M, model.N, observations.size(), observations.data(), \
				   model.pi.data(), model.A, model.B);
		}
		auto end = std::chrono::steady_clock::now();
		auto duration_us = std::chrono::duration_cast
			<std::chrono::milliseconds>(end - begin).count() / n_runs;

		time_list.emplace_back(duration_us);
	}

	Benchmark bench(time_list, impl_tag, "time (ms)", base_model.N, base_model.M, observations.size());
	bench.BeautyPrint(xout);
}


void perf_test_chrono(const std::string& impl_tag, compute_func baum_welch,
		int M, int N, int n_runs, int n_iter, std::ostream& xout) {

	std::vector<int> observations = uniform_emission_sample(SEQ_LEN, N);
	std::vector<double> time_list;

	HMM base_model(M, N);
	for (auto i = 0; i < n_iter; i++) {
		HMM model(base_model);

		auto begin = std::chrono::steady_clock::now();
		for (auto j = 0; j < n_runs; j++) {
			baum_welch(model.M, model.N, observations.size(), observations.data(), \
				   model.pi.data(), model.A, model.B);
		}
		auto end = std::chrono::steady_clock::now();
		auto duration_us = std::chrono::duration_cast
			<std::chrono::milliseconds>(end - begin).count() / n_runs;

		time_list.emplace_back(duration_us);
	}



	Benchmark bench(time_list, impl_tag, "time (ms)", base_model.N, base_model.M, observations.size());
	bench.BeautyPrint(xout);

}


void perf_test_chrono(const std::string& impl_tag, compute_func2 baum_welch,
		int M, int N, int n_runs, int n_iter, std::ostream& xout) {

	std::vector<int> observations = uniform_emission_sample(SEQ_LEN, N);
	std::vector<double> time_list;

	HMM base_model(M, N);
	for (auto i = 0; i < n_iter; i++) {
		HMM model(base_model);

		auto begin = std::chrono::steady_clock::now();
		for (auto j = 0; j < n_runs; j++) {
			baum_welch(model.transition, model.emission, model.pi, observations);
		}
		auto end = std::chrono::steady_clock::now();
		auto duration_us = std::chrono::duration_cast
			<std::chrono::milliseconds>(end - begin).count() / n_runs;

		time_list.emplace_back(duration_us);
	}

	Benchmark bench(time_list, impl_tag, "time (ms)", base_model.N, base_model.M, observations.size());
	bench.BeautyPrint(xout);

}
