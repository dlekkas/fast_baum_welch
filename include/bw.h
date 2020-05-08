#ifndef BW_H
#define BW_H

#include "hmm.h"

#define THRESHOLD 1e-4
#define MAX_ITERATIONS 5

void run_bw(int M, int N, int T, int* obs_sequence, double* pi, double** A, double** B,
		double** forward, double** backward, double** g, double*** chsi);

void run_bw_basic_opts(int M, int N, int T, int* obs_sequence, double* pi, double** A, double** B,
		double** forward, double** backward, double** g, double*** chsi);

void baum_welch(Matrix_v& transition, Matrix_v& emission, std::vector<double>& init_prob,
		const std::vector<int>& observations);

void baum_welch_opts(Matrix_v& transition, Matrix_v& emission, std::vector<double>& init_prob,
		const std::vector<int>& observation);


#endif
