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

    double sum = 0.0, acc = 0.0;


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
            acc = 0.0;
            for (int j = 0; j < M; j++) {
                acc += forward[t-1][j] * A[j][i];
			}
            forward[t][i] = B[observation_seq[t]][i] * acc;
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
            sum = 0.0;
            for (int j = 0; j < M; j++) {
                sum += backward[t+1][j] * A[i][j] * B[observation_seq[t+1]][j];
			}
            backward[t][i] = sum*sc_factors[t];
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
			int t;
			for (t = 0; t < T-4; t+=4) {
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

            double occurrences = 0.0;
            for (int t = 0; t < T; t++) {
                if (observation_seq[t] == vk) {
                    occurrences += (forward[t][i] * backward[t][i]) / sc_factors[t];
				}
            }
			B[vk][i] = occurrences / sum;
        }
    }


}

void run_bw_opts_v2(int M, int N, int T, int* obs_sequence, double* pi, double** A, double** B,
        double** forward, double** backward, double** g, double*** chsi) {

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
