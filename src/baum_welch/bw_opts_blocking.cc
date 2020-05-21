/*
* Straightforward implementation of the Baum-Welch
* https://en.wikipedia.org/wiki/Baum–Welch_algorithm
*/

/*
* BW OPTS (MANOS) + BLOCKING ON A MATRIX + BLOCKING ON B MATRIX (DOES NOT IMPROVE)
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

//WATCH OUT INLINING
// -----> TOTAL (M*((N+2)*T + M*(2T+1)) muls, M*(N*(T+1) + M + 2T) divs, M*(T*(N+1) + (T-1)*(M+1)) adds)
static void update_and_check(double** forward, double** backward, int M, int N, int T,
		double* pi, double** A, double** B, int* observation_seq, double *sc_factors, const vector<vector<int>>& obs_dict) { //double** sum, int blocking_size, int blocking_size_2) {

    double sum[4][4];  // CHANGE PARAMETER
    int blocking_size = 4; // CHANGE PARAMETER

    int blocking_size_2 = 4; // CHANGE PARAMETER

    for (int t = 0; t < T; t++) {
		  sc_factors[t] = 1.0 / sc_factors[t];
	  }

	 // pi computation
	  for (int i = 0; i < M; i++) {
		  pi[i] = forward[0][i] * backward[0][i] * sc_factors[0];
	  }

    // Compute A
    for (int i = 0; i < M; i += blocking_size) {
  		for (int j = 0; j < M; j += blocking_size) {
        // zero the sum_block

        for (int k = 0; k < blocking_size; k++) {
          for (int l = 0; l < blocking_size; l++) {
            sum[k][l] = 0.0;
          }
        }
        // iterate until last block since T-1 position is not included in computation
        int t = 0;
  	    for (; t < T - blocking_size - 1; t += blocking_size) {  // CHECK FOR SCALAR REPLACEMENT AND UNROLLING AT FINAL VERSION
          // computation of block on sum matrix
          for (int i_ = i; i_ < i + blocking_size; i_++) {
            int k = i_%blocking_size;   // strength reduction
            for (int j_ = j; j_ < j + blocking_size; j_++) {
              int l = j_%blocking_size;
              for (int t_ = t; t_ < t + blocking_size; t_++) {
                sum[k][l] += (backward[t_+1][j_] * B[observation_seq[t_+1]][j_]) * forward[t_][i_];
              }
            }
          }
        }
        // final block
        for (int i_ = i; i_ < i + blocking_size; i_++) {
          int k = i_%blocking_size;
          for (int j_ = j; j_ < j + blocking_size; j_++) {
            int l = j_%blocking_size;
            for (int t_ = t; t_ < t + blocking_size - 1; t_++) {
              sum[k][l] += (backward[t_+1][j_] * B[observation_seq[t_+1]][j_]) * forward[t_][i_];
            }
          }
        }
        // Write back to A matrix - common factor
        for (int i_ = i; i_ < i + blocking_size; i_++) {
          for (int j_ = j; j_ < j + blocking_size; j_++) {
            A[i_][j_] *= sum[i_%blocking_size][j_%blocking_size];
          }
        }
  		}
    }

    // Divide by sum of gammas
    for (int i = 0; i < M; i++) {
      double acc = 0.0;
      for (int t = 0; t < T-1; t++) {
        acc += (forward[t][i] * backward[t][i]) * sc_factors[t];
      }
      acc = 1.0 / acc;
      for (int j = 0; j < M; j++) {
        A[i][j] = A[i][j] * acc; // common factor
      }
    }

    for (int i = 0; i < M; i += blocking_size_2) {
      for (int j = 0; j < N; j += blocking_size_2) {
        for (int i_ = i; i_ < i + blocking_size_2; i_++) {
          for (int j_ = j; j_ < j + blocking_size_2; j_++) {
            double occurrences = 0.0;
            for (const auto& t: obs_dict[j_]) {
                occurrences += (forward[t][i_] * backward[t][i_]) * sc_factors[t];
      			}
            B[j_][i_] = occurrences;
          }
        }
      }
    }

    for (int i = 0; i < M; i++) {
      double acc_2 = 0.0;
      for (int t = 0; t < T; t++) {
        acc_2 += (forward[t][i] * backward[t][i]) * sc_factors[t];
  		}
      acc_2 = 1.0 / acc_2;
      for (int j = 0; j < N; j++) {
        B[j][i] = B[j][i] * acc_2;
      }
    }
}


void BaumWelchCOptsBlocking::operator()() {

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
