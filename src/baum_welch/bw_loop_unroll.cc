/*
* Straightforward implementation of the Baum-Welch
* https://en.wikipedia.org/wiki/Baum–Welch_algorithm
*/
#include <iostream>
#include <assert.h>
#include "../include/baum_welch.h"

using namespace std;

inline void forward_backward(double** forward, double** backward, int M, int N, int T,
		double* pi, double** A, double** B, int* observation_seq, double *sc_factors) {

    double sum = 0.0;

	// ops = 2*M , mem = 4*M
    for (int i = 0; i < M; i++) {
        forward[0][i] = pi[i] * B[observation_seq[0]][i];
        sum += forward[0][i];
    }
    sc_factors[0] = 1.0 / sum;

	// ops = M, mem = 2*M
    for (int i = 0; i < M; i++) {
        forward[0][i] *= sum;
    }

	// ops = 2*T*M^2 , mem = 2*T*M^2
    for (int t = 1; t < T; t++) {

        sum = 0.0;
	    for (int i = 0; i < M; i++) {
            double acc1 = 0.0, acc2 = 0.0, acc3 = 0.0, acc4 = 0.0;
			int j = 0;
            for (; j < M-3; j+=4) {
                acc1 += forward[t-1][j] * A[j][i];
                acc2 += forward[t-1][j+1] * A[j+1][i];
                acc3 += forward[t-1][j+2] * A[j+2][i];
                acc4 += forward[t-1][j+3] * A[j+3][i];
			}
            for (; j < M; j++) {
				acc1 += forward[t-1][j] * A[j][i];
			}
            forward[t][i] = B[observation_seq[t]][i] * (acc1 + acc2 + acc3 + acc4);
            sum += forward[t][i];
        }


        sc_factors[t] = 1.0 / sum;
        for (int i = 0; i < M; i++) {
            forward[t][i] = forward[t][i] / sum;
		}
    }

    for (int i = 0; i < M; i++) {
        backward[T-1][i] = sc_factors[T-1];
	}



	// ops = 3*T*M^2, mem = 3*T*M^2
    for (int t = T-2; t >= 0; t--) {
		int i = 0;
        for (; i < M-3; i+=4) {
            double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
			for (int j = 0; j < M; j+=1) {
				double tmp1 = backward[t+1][j] * B[observation_seq[t+1]][j];
				sum1 += A[i][j]   * tmp1;
				sum2 += A[i+1][j] * tmp1;
				sum3 += A[i+2][j] * tmp1;
				sum4 += A[i+3][j] * tmp1;
			}
            backward[t][i]   = sum1 * sc_factors[t];
			backward[t][i+1] = sum2 * sc_factors[t];
			backward[t][i+2] = sum3 * sc_factors[t];
			backward[t][i+3] = sum4 * sc_factors[t];
        }
		for (; i < M; i++) {
			double sum1 = 0.0;
			for (int j = 0; j < M; j++) {
				sum1 += A[i][j] * backward[t+1][j] * B[observation_seq[t+1]][j];
			}
			backward[t][i] = sum1 * sc_factors[t];
		}
    }

}


inline void update_and_check(double** forward, double** backward, int M, int N, int T,
		double* pi, double** A, double** B, int* observation_seq, double *sc_factors,
		const vector<vector<int>>& obs_dict) {

	for (int t = 0; t < T; t++) {
		sc_factors[t] = 1.0 / sc_factors[t];
	}

	for (int z = 0; z < M; z++) {
		pi[z] = forward[0][z] * backward[0][z] * sc_factors[0];
	}



	// ops = 3*T*M^2 , mem = 4*T*M^2
	int i = 0;
	for (; i < M-3; i+=4) {

        double acc1 = 0.0, acc2 = 0.0, acc3 = 0.0, acc4 = 0.0;
        for (int t = 0; t < T-1; t++) {
            acc1 += (forward[t][i] * backward[t][i]) * sc_factors[t];
            acc2 += (forward[t][i+1] * backward[t][i+1]) * sc_factors[t];
            acc3 += (forward[t][i+2] * backward[t][i+2]) * sc_factors[t];
            acc4 += (forward[t][i+3] * backward[t][i+3]) * sc_factors[t];
		}

		int j = 0;
		for (; j < M-1; j+=2) {
			double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
			double sum5 = 0.0, sum6 = 0.0, sum7 = 0.0, sum8 = 0.0;
			for (int t = 0; t < T-1; t++) {
				double c1 = backward[t+1][j]   * B[observation_seq[t+1]][j];
				double c2 = backward[t+1][j+1] * B[observation_seq[t+1]][j+1];
				sum1 += c1 * forward[t][i];   sum5 += c2 * forward[t][i];
				sum2 += c1 * forward[t][i+1]; sum6 += c2 * forward[t][i+1];
				sum3 += c1 * forward[t][i+2]; sum7 += c2 * forward[t][i+2];
				sum4 += c1 * forward[t][i+3]; sum8 += c2 * forward[t][i+3];
			}
			A[i][j]   *= sum1 / acc1;  A[i][j+1]   *= sum5 / acc1;
			A[i+1][j] *= sum2 / acc2;  A[i+1][j+1] *= sum6 / acc2;
			A[i+2][j] *= sum3 / acc3;  A[i+2][j+1] *= sum7 / acc3;
			A[i+3][j] *= sum4 / acc4;  A[i+3][j+1] *= sum8 / acc4;
		}

		for (; j < M; j++) {
			double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
			for (int t = 0; t < T-1; t++) {
				double c1 = backward[t+1][j] * B[observation_seq[t+1]][j];
				sum1 += c1 * forward[t][i];   sum2 += c1 * forward[t][i+1];
				sum3 += c1 * forward[t][i+2]; sum4 += c1 * forward[t][i+3];
			}
			A[i][j]   *= sum1 / acc1; A[i+1][j] *= sum2 / acc2;
			A[i+2][j] *= sum3 / acc3; A[i+3][j] *= sum4 / acc4;
		}
    }

	for (; i < M; i++) {
		double acc = 0.0;
		for (int t = 0; t < T-1; t++) {
			acc += (forward[t][i] * backward[t][i]) * sc_factors[t];
		}
		for (int j = 0; j < M; j++) {
			double sum = 0.0;
			for (int t = 0; t < T-1; t++) {
				sum += backward[t+1][j] * B[observation_seq[t+1]][j];
			}
			A[i][j] *= sum / acc;
		}
	}



	// ops <= 4*T*M*N, mem <= 3*T*M^2
    for (int k = 0; k < M; k++) {
        double sum1 = 0.0;
        for (int t = 0; t < T; t++) {
            sum1 += forward[t][k] * backward[t][k] * sc_factors[t];
		}

		int j = 0;
        for (; j < N-1; j+=2) {
            double acc1 = 0.0, acc2 = 0.0;
			for (const auto& t: obs_dict[j]) {
				acc1 += (forward[t][k] * backward[t][k]) * sc_factors[t];
			}
			for (const auto& t: obs_dict[j+1]) {
				acc2 += (forward[t][k] * backward[t][k]) * sc_factors[t];
			}
			B[j][k] = acc1 / sum1;
			B[j+1][k] = acc2 / sum1;
        }

		for (; j < N; j++) {
			double acc = 0.0;
			for (const auto& t: obs_dict[j]) {
				acc += (forward[t][k] * backward[t][k]) * sc_factors[t];
			}
			B[j][k] = acc / sum1;
		}
    }


}


void BaumWelchCLoopUnroll::operator()() {
	/* Significant optimization to get rid of conditional addition with
	 * extremely bad branch prediction patterns (highly unlikely to
	 * predict its outcome */
	vector<vector<int>> obs_dict(N);
	for (int t = 0; t < T; t++) {
		obs_dict[obs[t]].push_back(t);
	}


    double *sc_factors = (double *)malloc(T * sizeof(double));

	for (int i = 0; i < MAX_ITERATIONS; i++) {
        forward_backward(fwd, bwd, M, N, T, pi, A, B, obs, sc_factors);
        update_and_check(fwd, bwd, M, N, T, pi, A, B, obs, sc_factors, obs_dict);
    }

	#ifdef DEBUG
		cout << endl << "--------------------------  DEBUG START -------------------------- " << endl;
		cout << "Matrix A:" << endl;
		print_matrix(A, M, M);
		cout << endl << "Matrix B:" << endl;
		print_matrix(B, N, M);
		cout << endl << "Pi vector:" << endl;
		print_matrix(pi, M);
		cout << endl << "-------------------------- DEBUG END -------------------------- " << endl;
	#endif

}


inline void print_matrix(double** matrix, int n, int m) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
}

inline void print_vector(double* vec, int n) {
	for (int i = 0; i < n; i++) {
		cout << vec[i] << " ";
	}
	cout << endl;
}

