/*
* Straightforward implementation of the Baum-Welch
* https://en.wikipedia.org/wiki/Baum–Welch_algorithm
*/
#include <iostream>
#include <assert.h>
#include "../include/baum_welch.h"

using namespace std;

static void forward_backward(double** forward, double** backward, int M, int N, int T,
		double* pi, double** A, double** B, int* observation_seq, double *sc_factors) {

    double sum = 0.0;

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
	  // -> ops = 2*T*M^2 , mem = 2*T*M^2  ((T-1)*M*(M+1) muls, (T-1)*(M+1) divs, (T-1)*M*(M+1) adds)
    for (int t = 1; t < T; t++) {
      sum = 0.0;
			// -> M*(M+1) muls, M*(M+1) adds
	    for (int i = 0; i < M-3; i+=4) {
          double acc1 = 0.0, acc2 = 0.0, acc3 = 0.0, acc4 = 0.0;
					int j = 0;
          for (; j < M; j++) {
              acc1 += forward[t-1][j] * A[j][i];
              acc2 += forward[t-1][j] * A[j][i+1];
              acc3 += forward[t-1][j] * A[j][i+2];
              acc4 += forward[t-1][j] * A[j][i+3];
					}
          forward[t][i]   = B[observation_seq[t]][i]   * acc1;
          forward[t][i+1] = B[observation_seq[t]][i+1] * acc2;
          forward[t][i+2] = B[observation_seq[t]][i+2] * acc3;
          forward[t][i+3] = B[observation_seq[t]][i+3] * acc4;
          sum += (forward[t][i] + forward[t][i+1]);
          sum += (forward[t][i+2] + forward[t][i+3]);
      }
			// -> M+1 divs
      sc_factors[t] = 1.0 / sum;
      for (int i = 0; i < M; i++) {
          forward[t][i] = forward[t][i] * sc_factors[t];
			}
    }
		// -> not counted
    for (int i = 0; i < M; i++) {
        backward[T-1][i] = sc_factors[T-1];
		}
	// -> ops = 3*T*M^2, mem = 3*T*M^2   (EXACT FLOATS NOT COMPUTED)
  for (int t = T-2; t >= 0; t--) {
		int i = 0;
		// -> (M*(5M/4 +1) muls, M*M adds)
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


static void update_and_check(double** forward, double** backward, int M, int N, int T,
		double* pi, double** A, double** B, int* observation_seq, double *sc_factors,
		const vector<vector<int>>& obs_dict) {

  double sum[8][8];  // CHANGE PARAMETER
  int blocking_size = 8; // CHANGE PARAMETER

  int blocking_size_2 = 4; // CHANGE PARAMETER

	for (int t = 0; t < T; t++) {
		sc_factors[t] = 1.0 / sc_factors[t];
	}

	for (int z = 0; z < M; z++) {
		pi[z] = forward[0][z] * backward[0][z] * sc_factors[0];
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
      for (; t < T - blocking_size - 1; t += blocking_size) {


        // computation of block on sum matrix
        for (int i_ = i; i_ < i + blocking_size; i_ += 4) {
          int k = i_%blocking_size;   // strength reduction
          for (int j_ = j; j_ < j + blocking_size; j_ += 2) {
            int l = j_%blocking_size;
            double sum1 = 0.0, sum2 = 0.0, sum5 = 0.0, sum6 = 0.0;
            double sum3 = 0.0, sum4 = 0.0, sum7 = 0.0, sum8 = 0.0;

            //double sum9 = 0.0, sum10 = 0.0, sum11 = 0.0, sum12 = 0.0;
      			//double sum13 = 0.0, sum14 = 0.0, sum15 = 0.0, sum16 = 0.0;

            for (int t_ = t; t_ < t + blocking_size; t_++) {
              double c1 = backward[t_+1][j_]   * B[observation_seq[t_+1]][j_];
      				double c2 = backward[t_+1][j_+1] * B[observation_seq[t_+1]][j_+1];

             	//double c3 = backward[t_+1][j_+2]   * B[observation_seq[t_+1]][j_+2];
      				//double c4 = backward[t_+1][j_+3] * B[observation_seq[t_+1]][j_+3];

      				sum1 += c1 * forward[t_][i_];   sum5 += c2 * forward[t_][i_];
      				sum2 += c1 * forward[t_][i_+1]; sum6 += c2 * forward[t_][i_+1];
              sum3 += c1 * forward[t_][i_+2]; sum7 += c2 * forward[t_][i_+2];
      				sum4 += c1 * forward[t_][i_+3]; sum8 += c2 * forward[t_][i_+3];

              //sum9 += c3 * forward[t_][i_];   sum13 += c4 * forward[t_][i_];
      				//sum10 += c3 * forward[t_][i_+1]; sum14 += c4 * forward[t_][i_+1];
      				//sum11 += c3 * forward[t_][i_+2]; sum15 += c4 * forward[t_][i_+2];
      				//sum12 += c3 * forward[t_][i_+3]; sum16 += c4 * forward[t_][i_+3];
            }
            sum[k][l]   += sum1;  sum[k][l+1]   += sum5;
      			sum[k+1][l] += sum2;  sum[k+1][l+1] += sum6;
            sum[k+2][l] += sum3;  sum[k+2][l+1] += sum7;
      			sum[k+3][l] += sum4;  sum[k+3][l+1] += sum8;

  /*          sum[k][l+2]   += sum9;  sum[k][l+3]   += sum13;
      			sum[k+1][l+2] += sum10;  sum[k+1][l+3] += sum14;
      			sum[k+2][l+2] += sum11;  sum[k+2][l+3] += sum15;
      			sum[k+3][l+2] += sum12;  sum[k+3][l+3] += sum16; */
          }
        }


      }
      // final block
    /*  for (int i_ = i; i_ < i + blocking_size; i_++) {
        int k = i_%blocking_size;
        for (int j_ = j; j_ < j + blocking_size; j_++) {
          int l = j_%blocking_size;
          for (int t_ = t; t_ < t + blocking_size - 1; t_++) {
            sum[k][l] += (backward[t_+1][j_] * B[observation_seq[t_+1]][j_]) * forward[t_][i_];
          }
        }
      } */

			// final block
      for (int i_ = i; i_ < i + blocking_size; i_ += 4) {
        int k = i_%blocking_size;
        for (int j_ = j; j_ < j + blocking_size; j_ += 2) {
          int l = j_%blocking_size;
					double sum1 = 0.0, sum2 = 0.0, sum5 = 0.0, sum6 = 0.0;
					double sum3 = 0.0, sum4 = 0.0, sum7 = 0.0, sum8 = 0.0;
          for (int t_ = t; t_ < t + blocking_size - 1; t_++) {
						double c1 = backward[t_+1][j_]   * B[observation_seq[t_+1]][j_];
						double c2 = backward[t_+1][j_+1] * B[observation_seq[t_+1]][j_+1];

			/*      double c3 = backward[t_+1][j_+2]   * B[observation_seq[t_+1]][j_+2];
						double c4 = backward[t_+1][j_+3] * B[observation_seq[t_+1]][j_+3];*/

						sum1 += c1 * forward[t_][i_];   sum5 += c2 * forward[t_][i_];
						sum2 += c1 * forward[t_][i_+1]; sum6 += c2 * forward[t_][i_+1];
						sum3 += c1 * forward[t_][i_+2]; sum7 += c2 * forward[t_][i_+2];
						sum4 += c1 * forward[t_][i_+3]; sum8 += c2 * forward[t_][i_+3];
          }
					sum[k][l]   += sum1;  sum[k][l+1]   += sum5;
					sum[k+1][l] += sum2;  sum[k+1][l+1] += sum6;
					sum[k+2][l] += sum3;  sum[k+2][l+1] += sum7;
					sum[k+3][l] += sum4;  sum[k+3][l+1] += sum8;
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
  for (int i = 0; i < M; i ++) {
    double acc = 0.0;
    for (int t = 0; t < T-1; t++) {
      acc += (forward[t][i] * backward[t][i]) * sc_factors[t];
    }
    acc = 1.0 / acc;
    for (int j = 0; j < M; j++) {
      A[i][j] = A[i][j] * acc; // common factor
    }
  }

	/*for (int i = 0; i < M; i += 4) {
    double acc1 = 0.0, acc2 = 0.0, acc3 = 0.0, acc4 = 0.0;
    for (int t = 0; t < T-1; t++) {
			acc1 += (forward[t][i] * backward[t][i]) * sc_factors[t];
			acc2 += (forward[t][i+1] * backward[t][i+1]) * sc_factors[t];
			acc3 += (forward[t][i+2] * backward[t][i+2]) * sc_factors[t];
			acc4 += (forward[t][i+3] * backward[t][i+3]) * sc_factors[t];
    }
    acc1 = 1.0 / acc1; acc2 = 1.0 / acc2; acc3 = 1.0 / acc3; acc4 = 1.0 / acc4;
    for (int j = 0; j < M; j += 2) {
      A[i][j] = A[i][j] * acc1;   A[i][j+1] = A[i][j+1] * acc1; // common factor
			A[i+1][j] = A[i+1][j] * acc2;   A[i+1][j+1] = A[i+1][j+1] * acc2;
			A[i+2][j] = A[i+2][j] * acc3;   A[i+2][j+1] = A[i+2][j+1] * acc3;
			A[i+3][j] = A[i+3][j] * acc4;   A[i+3][j+1] = A[i+3][j+1] * acc4;

    }
  } */


	// ops <= 4*T*M*N, mem <= 3*T*M^2
  for (int k = 0; k < M; k++) {
    double sum1 = 0.0;
    for (int t = 0; t < T; t++) {
      sum1 += forward[t][k] * backward[t][k] * sc_factors[t];
		}
    sum1 = 1.0 / sum1;
		int j = 0;
    for (; j < N-1; j+=2) {
      double acc1 = 0.0, acc2 = 0.0;
			for (const auto& t: obs_dict[j]) {
				acc1 += (forward[t][k] * backward[t][k]) * sc_factors[t];
			}
			for (const auto& t: obs_dict[j+1]) {
			  acc2 += (forward[t][k] * backward[t][k]) * sc_factors[t];
			}
			B[j][k] = acc1 * sum1;
			B[j+1][k] = acc2 * sum1;
    }
		for (; j < N; j++) {
			double acc = 0.0;
			for (const auto& t: obs_dict[j]) {
				acc += (forward[t][k] * backward[t][k]) * sc_factors[t];
			}
			B[j][k] = acc * sum1;
		}
  }
}


void BaumWelchCBlocking_Unroll::operator()() {
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
