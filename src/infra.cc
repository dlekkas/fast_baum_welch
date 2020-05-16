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
#include "../include/bw_baseline.h"



void perf_test_rdtscp(const std::string& impl_tag, BaumWelch& impl,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV, std::string out_file) {

	std::vector<int> observations = uniform_emission_sample(S, N);
	std::vector<double> cycles_list;

	init_tsc();

	HMM base_model(M, N);
	for (auto i = 0; i < n_iter; i++) {
		impl.Load(base_model, observations);

		uint64_t start = start_tsc();
		for (auto j = 0; j < n_runs; j++) {
			impl();
		}
		uint64_t end = stop_tsc();
		uint64_t cycles = (end - start) / (double) n_runs;

		cycles_list.emplace_back(cycles);
	}

	Benchmark bench(cycles_list, impl_tag, "cycles", N, M, observations.size(), MAX_ITERATIONS);

	if (to_CSV) {
		bench.CSVPrint(out_file);
	} else {
		bench.CompactPrint(xout);
	}

}


void perf_test_chrono(const std::string& impl_tag, BaumWelch& impl,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV, std::string out_file) {

	std::vector<int> observations = uniform_emission_sample(S, N);
	std::vector<double> time_list;

	HMM base_model(M, N);
	for (auto i = 0; i < n_iter; i++) {
		impl.Load(base_model, observations);

		auto begin = std::chrono::steady_clock::now();
		for (auto j = 0; j < n_runs; j++) {
			impl();
		}
		auto end = std::chrono::steady_clock::now();
		auto duration_us = std::chrono::duration_cast
			<std::chrono::milliseconds>(end - begin).count() / n_runs;

		time_list.emplace_back(duration_us);
	}

	Benchmark bench(time_list, impl_tag, "msec", base_model.N, base_model.M,
			observations.size(), MAX_ITERATIONS);

	if (to_CSV) {
		bench.CSVPrint(out_file);
	} else {
		bench.CompactPrint(xout);
	}
}



bool IsValidImpl(BaumWelch& impl) {
	BaumWelchBaseline base_impl;
	int n = 64, m = 64, o = 128;

	HMM base_model(m, n);

	std::vector<int> obs = uniform_emission_sample(o, n);

	impl.Load(base_model, obs);
	impl();
	HMM test_hmm = impl.GetHMM();

	base_impl.Load(base_model, obs);
	base_impl();
	HMM base_hmm = base_impl.GetHMM();

	return base_hmm.IsSimilar(test_hmm);
}





double** allocate_2d(int M, int T) {

	double** ar = (double**)malloc(M * sizeof(double*));
    for (int i=0; i<M; i++)
        ar[i] = (double*)calloc(T, sizeof(double));


	printf("Address allocated: %p\n", ar);


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
