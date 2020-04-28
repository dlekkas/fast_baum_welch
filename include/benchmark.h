#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <iostream>

#include "hmm.h"

using compute_func = void (*)(int, int, int, int*, double*, double**, double**);

using compute_func2 = void (*)(Matrix_v&, Matrix_v&, std::vector<double>&,
		const std::vector<int>&);

void perf_test_rdtscp(const std::string& tag, const std::string& init_file,
		const std::string& obs_file, compute_func baum_welch, int n_runs,
		int n_iterations, std::ostream& os);

void perf_test_rdtscp(const std::string& tag, compute_func baum_welch,
		int M, int N, int n_runs, int n_iter, std::ostream& xout);

void perf_test_rdtscp(const std::string& tag, compute_func2 baum_welch,
		int M, int N, int n_runs, int n_iter, std::ostream& xout);

void perf_test_chrono(const std::string& tag, const std::string& init_file,
		const std::string& obs_file, compute_func2 baum_welch, int n_runs,
		int n_iterations, std::ostream& os);

void perf_test_chrono(const std::string& tag, compute_func baum_welch,
		int M, int N, int n_runs, int n_iter, std::ostream& xout);

void perf_test_chrono(const std::string& tag, compute_func2 baum_welch,
		int M, int N, int n_runs, int n_iter, std::ostream& xout);


#endif
