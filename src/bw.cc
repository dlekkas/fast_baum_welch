/*  
* Straightforward implementation of the Baum-Welch 
* https://en.wikipedia.org/wiki/Baum–Welch_algorithm
*/
#include <assert.h>
#include "../include/bw.h"
#include <iostream>

using namespace std;
double** g;

int it=0;

void forward_backward(double** forward, double** backward, int M, int N, int T, double* pi, double** A, double** B, int* observation_seq) {

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
        backward[i][T-1] = 1.0*sc_factors[T-1];

    for (t=T-2; t>=0; t--) {
        for (i=0; i<M; i++) {
            double sum = 0.0;
            for (j=0; j<M; j++)
                sum += backward[j][t+1] * A[i][j] * B[j][observation_seq[t+1]]; 
            backward[i][t] = sum*sc_factors[t];
        }
    }

}

bool update_and_check(double** forward, double** backward, int M, int N, int T, double* pi, double** A, double** B, int* observation_seq) {

    int i, j, k, w, t, vk;
    bool converged = true;

    for (t=0; t<T; t++) {
        double sum = 0.0;
        for (j=0; j<M; j++)
            sum += forward[j][t] * backward[j][t];
        for (i=0; i<M; i++) 
            g[i][t] = (forward[i][t] * backward[i][t])/sum;
    }

    double chsi[M][M][T];
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


    double new_pi[M];
    double new_A[M][M];
    double new_B[M][N];

    // estimate new initial vector, transition and emission matrixes

    for (i=0; i<M; i++)
        new_pi[i] = g[i][0];

    for (i=0; i<M; i++) {
        double sum = 0.0;
        for (t=0; t<T-1; t++)
            sum += g[i][t];
        for (j=0; j<M; j++) {
            double sum2 = 0.0;
            for (t=0; t<T-1; t++)
                sum2 += chsi[i][j][t];
            new_A[i][j] = sum2/sum;
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
            new_B[i][vk] = occurrences/sum;
        }
    }

    // compute difference from the old values

    double diff;
    double max_pi_diff = 0.0;
    for (i=0; i<M; i++) {
        diff = std::abs(pi[i] - new_pi[i]);
        if (diff > max_pi_diff)
            max_pi_diff = diff;

        pi[i] = new_pi[i]; // can the compiler reorder this with the previous instruction?
    }

    double max_A_diff = 0.0;
    for (i=0; i<M; i++) {
        for (j=0; j<M; j++) {
            diff = std::abs(A[i][j] - new_A[i][j]);
            if (diff > max_A_diff)
                max_A_diff = diff;

            A[i][j] = new_A[i][j];
        }
    }

    double max_B_diff = 0.0;
    for (i=0; i<M; i++) {
        for (j=0; j<N; j++) {
            diff = std::abs(B[i][j] - new_B[i][j]);
            if (diff > max_B_diff)
                max_B_diff = diff;

            B[i][j] = new_B[i][j];
        }
    }

    converged = (max_pi_diff < THRESHOLD) && (max_A_diff < THRESHOLD) && (max_B_diff < THRESHOLD);

    return converged;
}

void run_bw(int M, int N, int T, int* obs_sequence, double* pi, double** A, double** B) {

    double** forward = (double**)malloc(M * sizeof(double*));
    for (int i=0; i<M; i++)
        forward[i] = (double*)calloc(T, sizeof(double));

    double** backward = (double**)malloc(M * sizeof(double*));
    for (int i=0; i<M; i++)
        backward[i] = (double*)calloc(T, sizeof(double));

    g = (double**)malloc(M * sizeof(double*));
    for (int i=0; i<M; i++)
        g[i] = (double*)calloc(T, sizeof(double)); 


    bool has_converged = false;
    int iterations = 0;
    while (!has_converged && (iterations < MAX_ITERATIONS)) {
        forward_backward(forward, backward, M, N, T, pi, A, B, obs_sequence);
        has_converged = update_and_check(forward, backward, M, N, T, pi, A, B, obs_sequence);
        iterations++;
        it++;
    }

    //cout.precision(4);
    //cout << "This is new A: " << endl;
    for (int i=0; i<M; i++) {
        for (int j=0; j<M; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;


    //cout << "This is new B: " << endl;
    for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++) {
            cout << B[i][j] << " ";
        }
        cout << endl;
    }

    cout << endl;

    //cout << "This is new pi: " << endl;
    for (int i=0; i<M; i++) 
        cout << pi[i] << " ";
    cout << endl;

    cout << endl;

    for (int t=0; t<T; t++) {
        cout << 1- (g[0][t] > 0.5) << " ";
    }
    cout << endl;

}
