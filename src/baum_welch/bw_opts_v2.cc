/*
* Straightforward implementation of the Baum-Welch
* https://en.wikipedia.org/wiki/Baum–Welch_algorithm
*/

/*
* INLINED FUNCTIONS + ELIMINATED GAMMA AND CHSI
*/
#include <iostream>
#include <assert.h>
#include "../include/baum_welch.h"

using namespace std;

// -----> TOTAL (T-1)M(3M+2) + 2M muls, (T-1)(M+1) + 1 divs, (T-1)M(2*M + 1) + M adds (SAME COUNT AS BASIC_OPTS)
inline void forward_backward(double** forward, double** backward, int M, int N, int T,
		double* pi, double** A, double** B, int* observation_seq, double *sc_factors) {

    double sum = 0.0, acc = 0.0;


	 // -> ops = 2*M , mem = 4*M  (M muls, M adds)
    for (int i = 0; i < M; i++) {
        forward[0][i] = pi[i] * B[observation_seq[0]][i];
        sum += forward[0][i];
    }

    sc_factors[0] = 1.0 / sum;

	// -> ops = M, mem = 2*M  (M muls)
    for (int i = 0; i < M; i++) {
        forward[0][i] *= sum;
    }

	// -> ops = 2*T*M^2 , mem = 2*T*M^2   ((T-1)*M*(M+1) muls, (T-1)(M+1) divs, (T-1)*M*(M+1) adds)
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

		/*
		for (int i = 0; i < M; i++) {
			sum += forward[t][i];
		}
		*/

		sc_factors[t] = 1.0 / sum;
		for (int i = 0; i < M; i++) {
			forward[t][i] = forward[t][i] / sum;
		}
    }


	// -> not counted
  	for (int i = 0; i < M; i++) {
        backward[T-1][i] = sc_factors[T-1];
	}

	// -> ops = 3*T*M^2, mem = 3*T*M^2   ((T-1)*M*(2M+1) muls, (T-1)*M*M adds)
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

// -----> TOTAL (M*((N+2)*T + M*(2T+1)) muls, M*(N*(T+1) + M + 2T) divs, M*(T*(N+1) + (T-1)*(M+1)) adds)
inline void update_and_check(double** forward, double** backward, int M, int N, int T,
		double* pi, double** A, double** B, int* observation_seq, double *sc_factors) {

	// -> ops = 3*T*M^2 , mem = 4*T*M^2   (M*(M*(2T-1) + T) muls,  M*(M+T) divs, (T-1)*M*(M+1) adds)
	for (int i = 0; i < M; i++) {   // Compute A by rows

		// -> T muls, T divs, T-1 adds
		pi[i] = (forward[0][i] * backward[0][i]) / sc_factors[0];
    double acc = pi[i];
    for (int t = 1; t < T-1; t++) {
    	acc += (forward[t][i] * backward[t][i]) / sc_factors[t];
		}

		// -> M*(2T-1) muls, M divs, (T-1)*M adds
		for (int j = 0; j < M; j++) {
			double sum = 0.0;
			for (int t = 0; t < T-1; t++) {
          sum += (backward[t+1][j] * B[observation_seq[t+1]][j]) * forward[t][i];
      }
			A[i][j] *= sum / acc; // WHY IS THAT?
		}

  }

	// -> ops = 3*T*M*N, mem = 3*T*M^2   (T*M*(N+1) muls, M*N*(T+1) + M*T divs, T*M*(N+1) adds)
    for (int i = 0; i < M; i++) {
        double sum = 0.0;
				// -> T muls, T divs, T sums
        for (int t = 0; t < T; t++) {
            sum += (forward[t][i] / sc_factors[t]) * backward[t][i];
				}

				// -> N*T muls, N*(T+1) divs, N*T adds
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


void BaumWelchCOptsV2::operator()() {

    double *sc_factors = (double *)malloc(T * sizeof(double));

	for (int i = 0; i < MAX_ITERATIONS; i++) {
        forward_backward(fwd, bwd, M, N, T, pi, A, B, obs, sc_factors);
        update_and_check(fwd, bwd, M, N, T, pi, A, B, obs, sc_factors);
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
