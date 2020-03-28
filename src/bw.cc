/*  
* Straightforward implementation of the Baum-Welch 
* https://en.wikipedia.org/wiki/Baum–Welch_algorithm
*/
#include <assert.h>
#include "../include/bw.h"
#include <iostream>

using namespace std;

void BW::forward_backward(double** forward, double** backward) {

    int i, j, t;

    for (i=0; i<hmm->M; i++)
        forward[i][0] = hmm->pi[i] * hmm->B[i][observation_seq[0]];

    for (t=1; t<T; t++) {
	    for (i=0; i<hmm->M; i++) {
            double sum = 0.0;
            for (j=0; j<hmm->M; j++)
                sum += forward[j][t-1] * hmm->A[j][i];
            forward[i][t] = hmm->B[i][observation_seq[t]] * sum;
        }
    }

    for (i=0; i<hmm->M; i++)
        backward[i][T-1] = 1.0;

    for (t=T-2; t>=0; t--) {
        for (i=0; i<hmm->M; i++) {
            double sum = 0.0;
            for (j=0; j<hmm->M; j++)
                sum += backward[j][t+1] * hmm->A[i][j] * hmm->B[j][observation_seq[t+1]]; 
            backward[i][t] = sum;
        }
    }
}

bool BW::update_and_check(double** forward, double** backward) {

    int i, j, k, w, t, vk;
    double gamma[hmm->M][T];
    bool converged = true;

    for (t=0; t<T; t++) {
        double sum = 0.0;
        for (j=0; j<hmm->M; j++)
            sum += forward[j][t] * backward[j][t];
        for (i=0; i<hmm->M; i++) 
            gamma[i][t] = (forward[i][t] * backward[i][t])/sum;
    }

    double chsi[hmm->M][hmm->M][T];
    for (t=0; t<T-1; t++) {
        double sum = 0.0;
        for (k=0; k<hmm->M; k++) {
            for (w=0; w<hmm->M; w++) 
                sum += forward[k][t] * hmm->A[k][w] * backward[w][t+1] * hmm->B[w][observation_seq[t+1]];
        }
        assert(sum != 0.0);
        for (i=0; i<hmm->M; i++) 
            for (j=0; j<hmm->M; j++) {
                chsi[i][j][t] = forward[i][t] * hmm->A[i][j] * backward[j][t+1] * hmm->B[j][observation_seq[t+1]];
                chsi[i][j][t] = chsi[i][j][t]/sum;
            }
    }

    double new_pi[hmm->M];
    double new_A[hmm->M][hmm->M];
    double new_B[hmm->M][hmm->N];

    // estimate new initial vector, transition and emission matrixes

    for (i=0; i<hmm->M; i++)
        new_pi[i] = gamma[i][0];

    for (i=0; i<hmm->M; i++) {
        double sum = 0.0;
        for (t=0; t<T-1; t++)
            sum += gamma[i][t];
        for (j=0; j<hmm->M; j++) {
            double sum2 = 0.0;
            for (t=0; t<T-1; t++)
                sum2 += chsi[i][j][t];
            new_A[i][j] = sum2/sum;
        }
    }

    for (i=0; i<hmm->M; i++) {
        double sum = 0.0;
        for (t=0; t<T; t++)
            sum += gamma[i][t];
        for (vk=0; vk<hmm->N; vk++) {
            double occurrences  = 0.0;
            for (t=0; t<T; t++) {
                if (observation_seq[t] == vk)
                    occurrences += gamma[i][t];
            }
            new_B[i][vk] = occurrences/sum;
        }
    }

    // compute difference from the old values

    double diff;
    double max_pi_diff = 0.0;
    for (i=0; i<hmm->M; i++) {
        diff = std::abs(hmm->pi[i] - new_pi[i]);
        if (diff > max_pi_diff)
            max_pi_diff = diff;

        hmm->pi[i] = new_pi[i]; // can the compiler reorder this with the previous instruction?
    }

    double max_A_diff = 0.0;
    for (i=0; i<hmm->M; i++) {
        for (j=0; j<hmm->M; j++) {
            diff = std::abs(hmm->A[i][j] - new_A[i][j]);
            if (diff > max_A_diff)
                max_A_diff = diff;

            hmm->A[i][j] = new_A[i][j];
        }
    }

    double max_B_diff = 0.0;
    for (i=0; i<hmm->M; i++) {
        for (j=0; j<hmm->N; j++) {
            diff = std::abs(hmm->B[i][j] - new_B[i][j]);
            if (diff > max_B_diff)
                max_B_diff = diff;

            hmm->B[i][j] = new_B[i][j];
        }
    }

    converged = (max_pi_diff < threshold) && (max_A_diff < threshold) && (max_B_diff < threshold);

    return converged;
}

void BW::run_bw() {

    double** forward = (double**)malloc(hmm->M * sizeof(double*));
    for (int i=0; i<hmm->M; i++)
        forward[i] = (double*)calloc(T, sizeof(double));

    double** backward = (double**)malloc(hmm->M * sizeof(double*));
    for (int i=0; i<hmm->M; i++)
        backward[i] = (double*)calloc(T, sizeof(double));


    cout << "This is old pi: " << endl;
    for (int i=0; i<hmm->M; i++) 
        cout << hmm->pi[i] << " ";
    cout << endl;


    cout << "This is old A: " << endl;
    for (int i=0; i<hmm->M; i++) {
        for (int j=0; j<hmm->M; j++) {
            cout << hmm->A[i][j] << " ";
        }
        cout << endl;
    }

    cout << "This is old B: " << endl;
    for (int i=0; i<hmm->M; i++) {
        for (int j=0; j<hmm->N; j++) {
            cout << hmm->B[i][j] << " ";
        }
        cout << endl;
    }    

    bool has_converged = false;
    int iterations = 0;
    while (!has_converged && (iterations < max_iterations)) {
        forward_backward(forward, backward);
        has_converged = update_and_check(forward, backward);
        iterations++;
    }

    cout << "-----------------------------------------" << endl;

    cout << "This is new pi: " << endl;
    for (int i=0; i<hmm->M; i++) 
        cout << hmm->pi[i] << " ";
    cout << endl;


    cout << "This is new A: " << endl;
    for (int i=0; i<hmm->M; i++) {
        for (int j=0; j<hmm->M; j++) {
            cout << hmm->A[i][j] << " ";
        }
        cout << endl;
    }

    cout << "This is new B: " << endl;
    for (int i=0; i<hmm->M; i++) {
        for (int j=0; j<hmm->N; j++) {
            cout << hmm->B[i][j] << " ";
        }
        cout << endl;
    }



}
