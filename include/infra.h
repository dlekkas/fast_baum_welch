#ifndef INFRA_H
#define INFRA_H

#include <iostream>

#include "hmm.h"

using compute_func = void (*)(int, int, int, int*, double*, double**, double**,
		double**, double**, double**, double***);

using compute_func2 = void (*)(Matrix_v&, Matrix_v&, std::vector<double>&,
		const std::vector<int>&);

void perf_test_rdtscp(const std::string& tag, compute_func baum_welch,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV = false);

void perf_test_rdtscp(const std::string& tag, compute_func2 baum_welch,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV = false);

void perf_test_chrono(const std::string& tag, compute_func baum_welch,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV = false);

void perf_test_chrono(const std::string& tag, compute_func2 baum_welch,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV = false);


bool IsValidImpl(compute_func impl, compute_func base);

bool IsValidImpl(compute_func impl);



/* helper functions to preallocate/free memory needed for baum-welch */

double** allocate_2d(int M, int T);
double*** allocate_3d(int M, int K, int T);

double** free_2d(double** ar, int M, int T);
double*** free_3d(double*** ar, int M, int K, int T);

#endif
