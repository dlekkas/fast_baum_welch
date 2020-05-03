/*
* Straightforward implementation of the Baum-Welch
* https://en.wikipedia.org/wiki/Baum–Welch_algorithm
*/
#include <iostream>
#include <assert.h>
#include <math.h>
#include "../include/bw.h"

static int it=0;

void forward_backward(double** forward, double** backward, int M, int N, int T,
		double* pi, double** A, double** B, int* observation_seq) {


    int i, j, t;
    double sc_factors[T];

    double sum_i0 = 0.0;
    for (i=0; i<M; i++) {
        forward[i][0] = pi[i] * B[i][observation_seq[0]];
        sum_i0 += forward[i][0];
    }
    sc_factors[0] = 1.0/sum_i0;

    for (i=0; i<M; i++) {
        forward[i][0] = forward[i][0]*sc_factors[0];
    }

    for (t=1; t<T; t++) {
        double sum_i=0.0;
	    for (i=0; i<M; i++) {
            double sum = 0.0;
            for (j=0; j<M; j++)
                sum += forward[j][t-1] * A[j][i];
            forward[i][t] = B[i][observation_seq[t]] * sum;
            sum_i += forward[i][t];
        }
        sc_factors[t] = 1.0/sum_i;
        for (i=0; i<M; i++)
            forward[i][t] = forward[i][t]*sc_factors[t];
    }

    for (i=0; i<M; i++)
        backward[i][T-1] = sc_factors[T-1];

    for (t=T-2; t>=0; t--) {
        for (i=0; i<M; i++) {
            double sum = 0.0;
            for (j=0; j<M; j++)
                sum += backward[j][t+1] * A[i][j] * B[j][observation_seq[t+1]];
            backward[i][t] = sum*sc_factors[t];
        }
    }

}


bool update_and_check(double** forward, double** backward, int M, int N, int T,
		double* pi, double** A, double** B, int* observation_seq, double** g, double*** chsi) {

    int i, j, k, w, t, vk;
    bool converged = true;

    for (t=0; t<T; t++) {
        double sum = 0.0;
        for (j=0; j<M; j++)
            sum += forward[j][t] * backward[j][t];

        for (i=0; i<M; i++)
            g[i][t] = (forward[i][t] * backward[i][t])/sum;
    }

    for (t=0; t<T-1; t++) {
        double sum = 0.0;
        for (k=0; k<M; k++) {
            for (w=0; w<M; w++)
                sum += forward[k][t] * A[k][w] * backward[w][t+1] * B[w][observation_seq[t+1]];
        }
        assert(sum != 0.0);
        for (i=0; i<M; i++)
            for (j=0; j<M; j++) {
                chsi[i][j][t] = forward[i][t] * A[i][j] * backward[j][t+1] * B[j][observation_seq[t+1]];
                chsi[i][j][t] = chsi[i][j][t]/sum;
            }
    }


    double new_pi;
    double new_A;
    double new_B;

    // estimate new initial vector, transition and emission matrixes
    double diff;
    for (i=0; i<M; i++) {
        new_pi = g[i][0];
        diff = fabs(pi[i] - new_pi); // if we use mse the number of flops will change
        if (diff > THRESHOLD) {
            converged = false;
        }
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

            diff = fabs(A[i][j] - new_A);
            if (diff > THRESHOLD) {
                converged = false;
            }
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

            diff = fabs(B[i][vk] - new_B);
            if (diff > THRESHOLD) {
                converged = false;
            }

            B[i][vk] = new_B;
        }
    }

    return converged;
}

void run_bw(int M, int N, int T, int* obs_sequence, double* pi, double** A, double** B,
        double** forward, double** backward, double** g, double*** chsi) {

    int iterations = 0;
    while (iterations < MAX_ITERATIONS) { // SET TO MAX ITER
        forward_backward(forward, backward, M, N, T, pi, A, B, obs_sequence);
        update_and_check(forward, backward, M, N, T, pi, A, B, obs_sequence, g, chsi);
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
