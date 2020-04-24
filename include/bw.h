#ifndef BW_H
#define BW_H

#include "hmm.h"

#define THRESHOLD 1e-4
#define MAX_ITERATIONS 5

void forward_backward(double** forward, double** backward, int M, int N,
		int T, double* pi, double** A, double** B, int* observation_seq);

bool update_and_check(double** forward, double** backward, int M, int N,
		int T, double* pi, double** A, double** B, int* observation_seq);

void run_bw(int M, int N, int T, int* obs_sequence, double* pi, double** A, double** B);


#endif
