/*
* Straightforward implementation of the Baum-Welch
* https://en.wikipedia.org/wiki/Baum–Welch_algorithm
*/
#include <iostream>
#include <assert.h>
#include "../include/bw.h"

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
        for (int i = 0; i < M; i++) {
            double sum1 = 0.0, sum2 = 0.0;
			int j = 0;
            for (; j < M-1; j+=2) {
                sum1 += backward[t+1][j] * A[i][j] * B[observation_seq[t+1]][j];
                sum2 += backward[t+1][j+1] * A[i][j+1] * B[observation_seq[t+1]][j+1];
			}
			for (; j < M; j++) {
				sum1 += backward[t+1][j] * A[i][j] * B[observation_seq[t+1]][j];
			}
            backward[t][i] = (sum1+sum2)*sc_factors[t];
        }
    }

}


inline void update_and_check(double** forward, double** backward, int M, int N, int T,
		double* pi, double** A, double** B, int* observation_seq, double *sc_factors) {

	// ops = 3*T*M^2 , mem = 4*T*M^2
	for (int i = 0; i < M; i++) {

		pi[i] = (forward[0][i] * backward[0][i]) / sc_factors[0];

        double acc = pi[i];
        for (int t = 1; t < T-1; t++) {
            acc += (forward[t][i] * backward[t][i]) / sc_factors[t];
		}

		for (int j = 0; j < M; j++) {
			double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
			int t = 0;
			for (; t < T-4; t+=4) {
                sum1 += (backward[t+1][j] * B[observation_seq[t+1]][j]) * forward[t][i];
                sum2 += (backward[t+2][j] * B[observation_seq[t+2]][j]) * forward[t+1][i];
                sum3 += (backward[t+3][j] * B[observation_seq[t+3]][j]) * forward[t+2][i];
                sum4 += (backward[t+4][j] * B[observation_seq[t+4]][j]) * forward[t+3][i];
            }
			for (; t < T-1; t++) {
				sum1 += (backward[t+1][j] * B[observation_seq[t+1]][j]) * forward[t][i];
			}

			A[i][j] *= (sum1 + sum2 + sum3 + sum4) / acc;
		}

    }


	// ops = 3*T*M*N, mem = 3*T*M^2
    for (int i = 0; i < M; i++) {

        double sum = 0.0;
        for (int t = 0; t < T; t++) {
            sum += (forward[t][i] / sc_factors[t]) * backward[t][i];
		}

        for (int vk = 0; vk < N; vk++) {

            double acc1 = 0.0, acc2 = 0.0;
			int t = 0;
            for (; t < T-3; t+=4) {
                if (observation_seq[t] == vk) {
                    acc1 += (forward[t][i] * backward[t][i]) / sc_factors[t];
				}
                if (observation_seq[t+1] == vk) {
                    acc2 += (forward[t+1][i] * backward[t+1][i]) / sc_factors[t+1];
				}
                if (observation_seq[t+2] == vk) {
                    acc1 += (forward[t+2][i] * backward[t+2][i]) / sc_factors[t+2];
				}
                if (observation_seq[t+3] == vk) {
                    acc2 += (forward[t+3][i] * backward[t+3][i]) / sc_factors[t+3];
				}
            }
			for (; t < T; t++) {
				if (observation_seq[t] == vk) {
					acc1 += (forward[t][i] * backward[t][i]) / sc_factors[t];
				}
			}
			B[vk][i] = (acc1 + acc2) / sum;
        }
    }


}

void bw_loop_unroll(int M, int N, int T, int* obs_sequence, double* pi, double** A, double** B,
        double** forward, double** backward) {

    double *sc_factors = (double *)malloc(T * sizeof(double));

	for (int i = 0; i < MAX_ITERATIONS; i++) {
        forward_backward(forward, backward, M, N, T, pi, A, B, obs_sequence, sc_factors);
        update_and_check(forward, backward, M, N, T, pi, A, B, obs_sequence, sc_factors);
    }


	#ifdef DEBUG
		using namespace std;

		cout << endl << "--------------------------  DEBUG START -------------------------- " << endl;
		cout << "Matrix A:" << endl;
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				cout << A[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl << "Matrix B:" << endl;
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < N; j++) {
				cout << B[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl << "Pi vector:" << endl;
		for (int i = 0; i < M; i++) {
			cout << pi[i] << " ";
		}
		cout << endl << endl;
		cout << endl << "-------------------------- DEBUG END -------------------------- " << endl;
	#endif

}
