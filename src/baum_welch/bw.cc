/*
* Straightforward implementation of the Baum-Welch
* https://en.wikipedia.org/wiki/Baum–Welch_algorithm
*/
#include <iostream>
#include <assert.h>
#include <math.h>
#include "../include/baum_welch.h"

static int it=0;

// -----> TOTAL 3(T-1)M(M+1) + 2M muls, T divs, (T-1)M(2*M + 1) + M adds
void forward_backward(double** forward, double** backward, int M, int N, int T,
		double* pi, double** A, double** B, int* observation_seq) {


    int i, j, t;
    double sc_factors[T];

    double sum_i0 = 0.0;
		// -> M muls, M adds
    for (i = 0; i < M; i++) {
        forward[0][i] = pi[i] * B[observation_seq[0]][i];
        sum_i0 += forward[0][i];
    }
		// -> 1 div
    sc_factors[0] = 1.0/sum_i0;

		// -> M muls
    for (i = 0; i < M; i++) {
        forward[0][i] = forward[0][i] * sc_factors[0];
    }
		// -> (T-1)*M*(M+2) muls, T-1 divs, (T-1)*M*(M+1) adds
    for (t = 1; t < T; t++) {
      double sum_i = 0.0;
	    for (i = 0; i < M; i++) {
            double sum = 0.0;
            for (j = 0; j < M; j++)
                sum += forward[t-1][j] * A[j][i];
            forward[t][i] = B[observation_seq[t]][i] * sum;
            sum_i += forward[t][i];
      }
      sc_factors[t] = 1.0/sum_i;
      for (i = 0; i < M; i++)
          forward[t][i] = forward[t][i] * sc_factors[t];
    }

		// -> not counted
    for (i = 0; i < M; i++)
        backward[T-1][i] = sc_factors[T-1];

		// (T-1)*M*(2M+1) muls, (T-1)*M*M adds
    for (t = T-2; t >= 0; t--) {
        for (i = 0; i < M; i++) {
            double sum = 0.0;
            for (j = 0; j < M; j++)
                sum += backward[t+1][j] * A[i][j] * B[observation_seq[t+1]][j];
            backward[t][i] = sum * sc_factors[t];
        }
    }

}

// -----> TOTAL 2M((T-1)*3*M + T) muls, M*(T*(M+1) + N) divs, M*((T-1)*(2M+1) + T(N+2)) adds
bool update_and_check(double** forward, double** backward, int M, int N, int T,
		double* pi, double** A, double** B, int* observation_seq, double** g, double*** chsi) {

    int i, j, k, w, t, vk;
    bool converged = true;

		// -> T*2M muls, T*M divs, T*M adds
    for (t = 0; t < T; t++) {
        double sum = 0.0;
        for (j = 0; j < M; j++)
            sum += forward[t][j] * backward[t][j];

        for (i = 0; i < M; i++)
            g[i][t] = (forward[t][i] * backward[t][i])/sum;
    }

		// -> (T-1)*6*M*M muls, (T-1)*M*M divs, (T-1)*M*M adds
    for (t = 0; t < T-1; t++) {
        double sum = 0.0;
				// -> 3*M*M muls, M*M sums
        for (k = 0; k < M; k++) {
            for (w = 0; w < M; w++)
                sum += forward[t][k] * A[k][w] * backward[t+1][w] * B[observation_seq[t+1]][w];
        }
        assert(sum != 0.0);
				// -> 3*M*M muls, M*M divs
        for (i = 0; i < M; i++)
            for (j = 0; j < M; j++) {
                chsi[i][j][t] = forward[t][i] * A[i][j] * backward[t+1][j] * B[observation_seq[t+1]][j];
                chsi[i][j][t] = chsi[i][j][t]/sum;
            }
    }


    double new_pi;
    double new_A;
    double new_B;

		// -> not counted
    // estimate new initial vector, transition and emission matrixes
    for (i = 0; i < M; i++) {
        new_pi = g[i][0];
        pi[i] = new_pi;
    }

		// -> M*M divs, (T-1)*M*(M+1) adds
    for (i = 0; i < M; i++) {
        double sum = 0.0;
        for (t = 0; t < T-1; t++)
            sum += g[i][t];
        for (j = 0; j < M; j++) {
            double sum2 = 0.0;
            for (t = 0; t < T-1; t++)
                sum2 += chsi[i][j][t];
            new_A = sum2/sum;
            A[i][j] = new_A;
        }
    }

		// -> N*M divs, M*T*(N+1) adds
    for (i = 0; i < M; i++) {
        double sum = 0.0;
        for (t = 0; t < T; t++)
            sum += g[i][t];
        for (vk = 0; vk < N; vk++) {
            double occurrences  = 0.0;
            for (t = 0; t < T; t++) {
                if (observation_seq[t] == vk)
                    occurrences += g[i][t];
            }
            new_B = occurrences/sum;
            B[vk][i] = new_B;
        }
    }

    return converged;
}


void BaumWelchCBasic::operator()() {

    int iterations = 0;
    while (iterations < MAX_ITERATIONS) { // SET TO MAX ITER
        forward_backward(fwd, bwd, M, N, T, pi, A, B, obs);
        update_and_check(fwd, bwd, M, N, T, pi, A, B, obs, gamma, chsi);

        iterations++;
        it++;
    }

	#ifdef DEBUG
		using namespace std;

		cout << endl << "--------------------------  DEBUG START -------------------------- " << endl;

		cout << "Matrix A:" << endl;
		for (int i=0; i<M; i++) {
			for (int j=0; j<M; j++) {
				cout << A[i][j] << " ";
			}
			cout << endl;
		}

		cout << endl << "Matrix B:" << endl;
		for (int i=0; i<M; i++) {
			for (int j=0; j<N; j++) {
				cout << B[i][j] << " ";
			}
			cout << endl;
		}


		cout << endl << "Pi vector:" << endl;
		for (int i=0; i<M; i++)
			cout << pi[i] << " ";
		cout << endl << endl;

		for (int t=0; t<T; t++) {
			cout << 1- (g[0][t] > 0.5) << " ";
		}
		cout << endl << "-------------------------- DEBUG END -------------------------- " << endl;

	#endif

}
