#ifndef INFRA_H
#define INFRA_H

#include <iostream>

#include "hmm.h"
#include "bw.h"


using compute_func  = decltype(&run_bw);
using compute_func1 = decltype(&bw_loop_unroll);
using compute_func2 = decltype(&baum_welch);


void perf_test_rdtscp(const std::string& tag, compute_func baum_welch,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV = false,
		std::string out_file = "");

void perf_test_rdtscp(const std::string& tag, compute_func1 baum_welch,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV = false,
		std::string out_file = "");

void perf_test_rdtscp(const std::string& tag, compute_func2 baum_welch,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV = false,
		std::string out_file = "");

void perf_test_chrono(const std::string& tag, compute_func baum_welch,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV = false,
		std::string out_file = "");

void perf_test_chrono(const std::string& tag, compute_func1 baum_welch,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV = false,
		std::string out_file = "");

void perf_test_chrono(const std::string& tag, compute_func2 baum_welch,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV = false,
		std::string out_file = "");


bool IsValidImpl(compute_func impl);

bool IsValidImpl(compute_func1 impl);

bool IsValidImpl(compute_func2 impl);


/* helper functions to preallocate/free memory needed for baum-welch */

double** allocate_2d(int M, int T);
double*** allocate_3d(int M, int K, int T);

double** free_2d(double** ar, int M, int T);
double*** free_3d(double*** ar, int M, int K, int T);

#endif
