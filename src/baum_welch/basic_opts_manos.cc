/*
* Straightforward implementation of the Baum-Welch
* https://en.wikipedia.org/wiki/Baum–Welch_algorithm
*/

/*
* INLINED FUNCTIONS (DID NOT IMPROVE) + ELIMINATED GAMMA AND CHSI (LESS MEMORY - MORE COMPUTATIONS) + COMMON FACTOR IN A COMPUTATION + ELIMINATED CONDITION (OBS_DICT) + FEWER DIVS-MORE MULTS WHEN SCALING
*/
#include <iostream>
#include <assert.h>
#include "../include/baum_welch.h"

using namespace std;

// -----> TOTAL (T-1)M(3M+2) + 2M muls, (T-1)(M+1) + 1 divs, (T-1)M(2*M + 1) + M adds (SAME COUNT AS BASIC_OPTS)
static void forward_backward(double** forward, double** backward, int M, int N, int T,
		double* pi, double** A, double** B, int* observation_seq, double *sc_factors) {

    double sum = 0.0, acc = 0.0;

	 // -> ops = 2*M , mem = 4*M  (M muls, M adds)
    for (int i = 0; i < M; i++) {
    	forward[0][i] = pi[i] * B[observation_seq[0]][i];
      sum += forward[0][i];
    }
	// -> 1 div
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
      sc_factors[t] = 1.0 / sum;
      for (int i = 0; i < M; i++) {
        forward[t][i] = forward[t][i] * sc_factors[t];
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
          backward[t][i] = sum * sc_factors[t];
      }
    }
}

// -----> TOTAL (M*((N+2)*T + M*(2T+1)) muls, M*(N*(T+1) + M + 2T) divs, M*(T*(N+1) + (T-1)*(M+1)) adds)
static void update_and_check(double** forward, double** backward, int M, int N, int T,
		double* pi, double** A, double** B, int* observation_seq, double *sc_factors, const vector<vector<int>>& obs_dict) {

	for (int t = 0; t < T; t++) {
		sc_factors[t] = 1.0 / sc_factors[t];
	}

	// pi computation
	for (int i = 0; i < M; i++) {
		pi[i] = forward[0][i] * backward[0][i] * sc_factors[0];
	}

	// A computation
	// -> ops = 3*T*M^2 , mem = 4*T*M^2   (M*(M*(2T-1) + T) muls,  M*(M+T) divs, (T-1)*M*(M+1) adds)
	for (int i = 0; i < M; i++) {
		double acc = 0.0;
		for (int t = 0; t < T-1; t++) {
      acc += (forward[t][i] * backward[t][i]) * sc_factors[t];
    }
		acc = 1.0 / acc;
		for (int j = 0; j < M; j ++) {
			double sum = 0.0;
			for (int t = 0; t < T-1; t ++) {
          sum += (backward[t+1][j] * B[observation_seq[t+1]][j]) * forward[t][i];
      }
			A[i][j] *= sum * acc; // common factor
		}
  }

	// B computation
 	// -> ops = 3*T*M*N, mem = 3*T*M^2   (T*M*(N+1) muls, M*(N*(T+1) + T) divs, T*M*(N+1) adds)
  for (int i = 0; i < M; i++) {
    double acc_2 = 0.0;
    for (int t = 0; t < T; t++) {
        acc_2 += (forward[t][i] * backward[t][i]) * sc_factors[t];
		}
		acc_2 = 1.0 / acc_2;
    for (int j = 0; j < N; j++) {
        double occurrences = 0.0;
        for (const auto& t: obs_dict[j]) {
            occurrences += (forward[t][i] * backward[t][i]) * sc_factors[t];
				}
				B[j][i] = occurrences * acc_2;
    }
  }

	// Denominators computed at the end -> Same runtime
	/*for (int i = 0; i < M; i++) {
		double acc = 0.0;
		double acc_2;
		for (int t = 0; t < T-1; t++) {
			acc += (forward[t][i] * backward[t][i]) * sc_factors[t];
		}
		acc_2 = acc + (forward[T-1][i] * backward[T-1][i]) * sc_factors[T-1];
		acc = 1.0 / acc;

		acc_2 = 1.0 / acc_2;

		for (int j = 0; j < M; j++) {
			A[i][j] *= acc;
		}

		for (int j = 0; j < N; j++) {
			B[j][i] *= acc_2;
		}
	} */
}


void BaumWelchCOptsManos::operator()() {

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
