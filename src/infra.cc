#include <chrono>
#include <iostream>
#include <fstream>
#include <numeric>
#include <iterator>
#include <algorithm>

#include "../include/tsc_x86.h"
#include "../include/benchmark.h"
#include "../include/infra.h"


void perf_test_rdtscp(const std::string& impl_tag, compute_func baum_welch,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV) {

	std::vector<int> observations = uniform_emission_sample(S, N);
	std::vector<double> cycles_list;

	init_tsc();

	HMM base_model(M, N);
	for (auto i = 0; i < n_iter; i++) {

		HMM model(base_model);
		double** forward = allocate_2d(model.M, observations.size());
		double** backward = allocate_2d(model.M, observations.size());
		double** g = allocate_2d(model.M, observations.size());
		double*** chsi = allocate_3d(model.M, model.M, observations.size());

		uint64_t start = start_tsc();
		for (auto j = 0; j < n_runs; j++) {
			baum_welch(model.M, model.N, observations.size(), observations.data(), \
				   model.pi.data(), model.A, model.B, forward, backward, g, chsi);
		}
		uint64_t end = stop_tsc();
		uint64_t cycles = (end - start) / (double) n_runs;
		forward = free_2d(forward, model.M, observations.size());
		backward = free_2d(backward, model.M, observations.size());
		g = free_2d(g, model.M, observations.size());
		chsi = free_3d(chsi, model.M, model.M, observations.size());

		cycles_list.emplace_back(cycles);
	}


	Benchmark bench(cycles_list, impl_tag, "cycles", N, M, observations.size());

	if (to_CSV) {
		bench.CSVPrint("results_cycles.txt");
	} else {
		bench.CompactPrint(xout);
	}

}


void perf_test_rdtscp(const std::string& impl_tag, compute_func2 baum_welch,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV) {

	std::vector<int> observations = uniform_emission_sample(S, N);
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

	if (to_CSV) {
		bench.CSVPrint("results_cycles.txt");
	} else {
		bench.CompactPrint(xout);
	}


}




void perf_test_chrono(const std::string& impl_tag, compute_func baum_welch,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV) {

	std::vector<int> observations = uniform_emission_sample(S, N);
	std::vector<double> time_list;

	HMM base_model(M, N);
	for (auto i = 0; i < n_iter; i++) {
		HMM model(base_model);

		double** forward = allocate_2d(model.M, observations.size());
		double** backward = allocate_2d(model.M, observations.size());
		double** g = allocate_2d(model.M, observations.size());
		double*** chsi = allocate_3d(model.M, model.M, observations.size());

		auto begin = std::chrono::steady_clock::now();
		for (auto j = 0; j < n_runs; j++) {
			baum_welch(model.M, model.N, observations.size(), observations.data(), \
				   model.pi.data(), model.A, model.B, forward, backward, g, chsi);
		}
		auto end = std::chrono::steady_clock::now();
		auto duration_us = std::chrono::duration_cast
			<std::chrono::milliseconds>(end - begin).count() / n_runs;

		forward = free_2d(forward, model.M, observations.size());
		backward = free_2d(backward, model.M, observations.size());
		g = free_2d(g, model.M, observations.size());
		chsi = free_3d(chsi, model.M, model.M, observations.size());

		time_list.emplace_back(duration_us);
	}



	Benchmark bench(time_list, impl_tag, "msec", base_model.N, base_model.M,
			observations.size());

	if (to_CSV) {
		bench.CSVPrint("results_time.txt");
	} else {
		bench.CompactPrint(xout);
	}

}


void perf_test_chrono(const std::string& impl_tag, compute_func2 baum_welch,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV) {

	std::vector<int> observations = uniform_emission_sample(S, N);
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

	Benchmark bench(time_list, impl_tag, "msec", base_model.N, base_model.M,
			observations.size());

	if (to_CSV) {
		bench.CSVPrint("results_time.txt");
	} else {
		bench.CompactPrint(xout);
	}
}

double** allocate_2d(int M, int T) {

	double** ar = (double**)malloc(M * sizeof(double*));
    for (int i=0; i<M; i++)
        ar[i] = (double*)calloc(T, sizeof(double));
    return ar;
}

double*** allocate_3d(int M, int K, int T) {

	double*** ar = (double***)malloc(M * sizeof(double**));
	for (int i=0; i<M; i++) {
        ar[i] = (double**)malloc(K * sizeof(double*));
        for (int j=0; j<K; j++)
	        ar[i][j] = (double*)calloc(T, sizeof(double));
	}
	return ar;
}

double** free_2d(double** ar, int M, int T) {

    for (int i=0; i<M; i++)
        free(ar[i]);
    free(ar);
	return NULL;
}

double*** free_3d(double*** ar, int M, int K, int T) {

    for (int i=0; i<M; i++) {
		for (int j=0; j<K; j++)
			free(ar[i][j]);
		free(ar[i]);
	}
	free(ar);
	return NULL;
}
