/*
* Straightforward implementation of the Baum-Welch
* https://en.wikipedia.org/wiki/Baum–Welch_algorithm
*/
#include <iostream>
#include <assert.h>
#include "../include/baum_welch.h"

static int it=0;

static void forward_backward(double** forward, double** backward, int M, int N, int T,
		double* pi, double** A, double** B, int* observation_seq, double *sc_factors) {

    int i, j, t;

    double sum_i0 = 0.0;
    for (i=0; i<M; i++) {
        forward[0][i] = pi[i] * B[observation_seq[0]][i];
        sum_i0 += forward[0][i];
    }
    sc_factors[0] = 1.0/sum_i0;

    for (i=0; i<M; i++) {
        forward[0][i] = forward[0][i]*sc_factors[0];
    }

    for (t=1; t<T; t++) {
        double sum_i=0.0;
	    for (i=0; i<M; i++) {
            double sum = 0.0;
            for (j=0; j<M; j++)
                sum += forward[t-1][j] * A[j][i];
            forward[t][i] = B[observation_seq[t]][i] * sum;
            sum_i += forward[t][i];
        }
        sc_factors[t] = 1.0/sum_i;
        for (i=0; i<M; i++)
            forward[t][i] = forward[t][i]*sc_factors[t];
    }

    for (i=0; i<M; i++)
        backward[T-1][i] = sc_factors[T-1];

    for (t=T-2; t>=0; t--) {
        for (i=0; i<M; i++) {
            double sum = 0.0;
            for (j=0; j<M; j++)
                sum += backward[t+1][j] * A[i][j] * B[observation_seq[t+1]][j];
            backward[t][i] = sum*sc_factors[t];
        }
    }

}


static bool update_and_check(double** forward, double** backward, int M, int N, int T,
		double* pi, double** A, double** B, int* observation_seq, double *sc_factors, double** g, double*** chsi) {

    int i, j, t, vk;
    bool converged = true;

    // denominators of gamma[i][t] and chsi[i][j][t] are the same
    // and independent of i, j
    // double denom[T];

    for (t=0; t<T; t++) {
        for (i=0; i<M; i++)
            g[i][t] = (forward[t][i] * backward[t][i])/sc_factors[t];
    }

    //double den;
    for (t=0; t<T-1; t++) {
        // We proved by induction that the denominator of this chsi[i][j][t] is 1. We need to validate it and review the proof to be sure.
        // If this fails, we can merge the computation of gamma and chsi.
        for (i=0; i<M; i++)
            for (j=0; j<M; j++) {
                chsi[i][j][t] = forward[t][i] * A[i][j] * backward[t+1][j] * B[observation_seq[t+1]][j];
            }
    }


    double new_pi;
    double new_A;
    double new_B;

    // estimate new initial vector, transition and emission matrixes
    for (i=0; i<M; i++) {
        new_pi = g[i][0];
        pi[i] = new_pi;
    }

    for (i=0; i<M; i++) {
        double sum = 0.0;
        for (t=0; t<T-1; t++)
            sum += g[i][t];
        for (j=0; j<M; j++) {
            double sum2 = 0.0;
            for (t=0; t<T-1; t++)
                sum2 += chsi[i][j][t];
            new_A = sum2/sum;
            A[i][j] = new_A;

        }
    }

    for (i=0; i<M; i++) {
        double sum = 0.0;
        for (t=0; t<T; t++)
            sum += g[i][t];
        for (vk=0; vk<N; vk++) {
            double occurrences  = 0.0;
            for (t=0; t<T; t++) {
                if (observation_seq[t] == vk)
                    occurrences += g[i][t];
            }
            new_B = occurrences/sum;

            B[vk][i] = new_B;
        }
    }

    return converged;
}


void BaumWelchCBasicOpts::operator()() {

    int iterations = 0;
    double *sc_factors = (double *)malloc(T * sizeof(double));

    while (iterations < MAX_ITERATIONS) {
        forward_backward(fwd, bwd, M, N, T, pi, A, B, obs, sc_factors);
        update_and_check(fwd, bwd, M, N, T, pi, A, B, obs, sc_factors, gamma, chsi);
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
