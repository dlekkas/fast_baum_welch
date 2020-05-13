#include <chrono>
#include <iostream>
#include <fstream>
#include <numeric>
#include <iterator>
#include <algorithm>

#include "../include/tsc_x86.h"
#include "../include/benchmark.h"
#include "../include/infra.h"
#include "../include/bw.h"


void perf_test_rdtscp(const std::string& impl_tag, compute_func baum_welch,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV) {

	std::vector<int> observations = uniform_emission_sample(S, N);
	std::vector<double> cycles_list;

	init_tsc();

	HMM base_model(M, N);
	for (auto i = 0; i < n_iter; i++) {

		HMM model(base_model);
		double** forward = allocate_2d(observations.size(), model.M);
		double** backward = allocate_2d(observations.size(), model.M);
		double** g = allocate_2d(model.M, observations.size());
		double*** chsi = allocate_3d(model.M, model.M, observations.size());

		uint64_t start = start_tsc();
		for (auto j = 0; j < n_runs; j++) {
			baum_welch(model.M, model.N, observations.size(), observations.data(), \
				   model.pi.data(), model.A, model.B, forward, backward, g, chsi);
		}
		uint64_t end = stop_tsc();
		uint64_t cycles = (end - start) / (double) n_runs;
		forward = free_2d(forward, observations.size(), model.M);
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


void perf_test_rdtscp(const std::string& impl_tag, compute_func1 baum_welch,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV) {

	std::vector<int> observations = uniform_emission_sample(S, N);
	std::vector<double> cycles_list;

	init_tsc();

	HMM base_model(M, N);
	for (auto i = 0; i < n_iter; i++) {

		HMM model(base_model);
		double** forward = allocate_2d(observations.size(), model.M);
		double** backward = allocate_2d(observations.size(), model.M);

		uint64_t start = start_tsc();
		for (auto j = 0; j < n_runs; j++) {
			baum_welch(model.M, model.N, observations.size(), observations.data(), \
				   model.pi.data(), model.A, model.B, forward, backward);
		}
		uint64_t end = stop_tsc();
		uint64_t cycles = (end - start) / (double) n_runs;
		forward = free_2d(forward, observations.size(), model.M);
		backward = free_2d(backward, model.M, observations.size());

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

		double** forward = allocate_2d(observations.size(), model.M);
		double** backward = allocate_2d(observations.size(), model.M);
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

		forward = free_2d(forward, observations.size(), model.M);
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


void perf_test_chrono(const std::string& impl_tag, compute_func1 baum_welch,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV) {

	std::vector<int> observations = uniform_emission_sample(S, N);
	std::vector<double> time_list;

	HMM base_model(M, N);
	for (auto i = 0; i < n_iter; i++) {
		HMM model(base_model);

		double** forward = allocate_2d(observations.size(), model.M);
		double** backward = allocate_2d(observations.size(), model.M);

		auto begin = std::chrono::steady_clock::now();
		for (auto j = 0; j < n_runs; j++) {
			baum_welch(model.M, model.N, observations.size(), observations.data(), \
				   model.pi.data(), model.A, model.B, forward, backward);
		}
		auto end = std::chrono::steady_clock::now();
		auto duration_us = std::chrono::duration_cast
			<std::chrono::milliseconds>(end - begin).count() / n_runs;

		forward = free_2d(forward, observations.size(), model.M);
		backward = free_2d(backward, model.M, observations.size());

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


bool IsValidImpl(compute_func impl) {
	int n = 64, m = 64, o = 128;
	HMM base_model(m, n), test_model(base_model);

	std::vector<int> obs = uniform_emission_sample(o, n);
	baum_welch(base_model.transition, base_model.emission, base_model.pi, obs);
	for (auto i = 0; i < m; i++) {
		std::copy(base_model.transition[i].begin(), base_model.transition[i].end(), base_model.A[i]);
	}
	for (auto i = 0; i < n; i++) {
		std::copy(base_model.emission[i].begin(), base_model.emission[i].end(), base_model.B[i]);
	}

	double** fwd = allocate_2d(o, m);
	double** bwd = allocate_2d(o, m);
	double** g = allocate_2d(m, o);
	double*** chsi = allocate_3d(m, m, o);
	impl(m, n, o, obs.data(), test_model.pi.data(), test_model.A, test_model.B, fwd, bwd, g, chsi);

	return test_model.IsSimilar(base_model);
}


bool IsValidImpl(compute_func1 impl) {
	int n = 64, m = 64, o = 128;
	HMM base_model(m, n), test_model(base_model);

	std::vector<int> obs = uniform_emission_sample(o, n);
	baum_welch(base_model.transition, base_model.emission, base_model.pi, obs);
	for (auto i = 0; i < m; i++) {
		std::copy(base_model.transition[i].begin(), base_model.transition[i].end(), base_model.A[i]);
	}
	for (auto i = 0; i < n; i++) {
		std::copy(base_model.emission[i].begin(), base_model.emission[i].end(), base_model.B[i]);
	}

	double** fwd = allocate_2d(o, m);
	double** bwd = allocate_2d(o, m);
	impl(m, n, o, obs.data(), test_model.pi.data(), test_model.A, test_model.B, fwd, bwd);

	return test_model.IsSimilar(base_model);
}


bool IsValidImpl(compute_func2 impl) {
	int n = 64, m = 64, o = 128;

	HMM base_model(m, n);
	HMM test_model(base_model);

	std::vector<int> obs = uniform_emission_sample(o, n);

	baum_welch(base_model.transition, base_model.emission, base_model.pi, obs);
	for (auto i = 0; i < m; i++) {
		std::copy(base_model.transition[i].begin(), base_model.transition[i].end(), base_model.A[i]);
	}

	for (auto i = 0; i < n; i++) {
		std::copy(base_model.emission[i].begin(), base_model.emission[i].end(), base_model.B[i]);
	}


	impl(test_model.transition, test_model.emission, test_model.pi, obs);
	for (auto i = 0; i < m; i++) {
		std::copy(test_model.transition[i].begin(), test_model.transition[i].end(), test_model.A[i]);
	}

	for (auto i = 0; i < n; i++) {
		std::copy(test_model.emission[i].begin(), test_model.emission[i].end(), test_model.B[i]);
	}

	return test_model.IsSimilar(base_model);
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
