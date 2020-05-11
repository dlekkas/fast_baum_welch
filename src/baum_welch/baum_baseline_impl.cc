#include <iostream>
#include <vector>

#include "../include/bw.h"

#define MAX_ITER 5

using namespace std;

using Matrix_v = vector<vector<double>>;

static Matrix_v forward(const Matrix_v&, const Matrix_v&, const vector<double>&,
		const vector<int>&, vector<double>&);

static Matrix_v backward(const Matrix_v&,  const Matrix_v&, const vector<double>&,
		const vector<int>&, const vector<double>&);

static void maximization_step(Matrix_v&, Matrix_v&, vector<double>&,
		const vector<int>&, const Matrix_v&, const vector<Matrix_v>&);

static void expectation_step(const Matrix_v&, const Matrix_v&, const Matrix_v&,
		const Matrix_v&, const vector<int>&, Matrix_v&, vector<Matrix_v>&);





void baum_welch(Matrix_v& transition, Matrix_v& emission, vector<double>& init_prob,
		const vector<int>& observation) {

	int T = observation.size();
	int n_states = transition.size();



	// probability of being in state `i` at time `t` and state `j` at time `t+1`,
	// given the observation sequence and HMM model
	vector<Matrix_v> ksi(T, Matrix_v(n_states, vector<double>(n_states)));

	// probability of being in state `i` at time `t`
	Matrix_v gamma(n_states, vector<double>(T));

	// scaling factors required to avoid summation to zero due to limited precision
	vector<double> scale_c(T);

	for (int i = 0; i < MAX_ITERATIONS; i++) {

		Matrix_v fwd = forward(transition, emission, init_prob, observation, scale_c);
		Matrix_v bwd = backward(transition, emission, init_prob, observation, scale_c);

		// EXPECTATION step
		expectation_step(transition, emission, fwd, bwd, observation, gamma, ksi);

		// MAXIMIZATION step
		maximization_step(transition, emission, init_prob, observation, gamma, ksi);
	}


	#ifdef DEBUG
		cout << endl << "Transition matrix" << endl;
		for (auto row: transition) {
			for (auto x: row) { cout << x << " "; }
			cout << endl;
		}

		cout << endl << "Emission matrix" << endl;
		for (auto row: emission) {
			for (auto x: row) { cout << x << " "; }
			cout << endl;
		}

		cout << endl << "Init prob vector" << endl;
		for (auto x: init_prob) { cout << x << " "; }
		cout << endl;
	#endif
}



/* calculate the expected state occupancy count and the expected
 * state transition count based on the HMM parameters of previous step */
void expectation_step(const Matrix_v& transition, const Matrix_v& emission,
		const Matrix_v& fwd, const Matrix_v& bwd, const vector<int>& observation,
		Matrix_v& gamma, vector<Matrix_v>& ksi) {

	int n_states = transition.size();
	int T = observation.size();

	for (int t = 0; t < T-1; t++) {
		// probability of observing sequence `observation` given the HMM model
		double obs_prob = 0.0;
		for (int j = 0; j < n_states; j++) {
			obs_prob += fwd[t][j] * bwd[t][j];
		}

		for (int i = 0; i < n_states; i++) {

			for (int j = 0; j < n_states; j++) {
				// probability of being in state `i` at time `t`, at state `j`
				// at time `t+1` and observing sequence `observation`
				double joint_prob = fwd[t][i] * transition[i][j] *
					emission[observation[t+1]][j] * bwd[t+1][j];

				// straightforward derivation from Bayes' theorem
				ksi[t][i][j] = joint_prob ;/// obs_prob;
			}
		}
	}


	for (int t = 0; t < T; t++) {
		// probability of observing sequence `observation` given the HMM model
		double obs_prob = 0.0;
		for (int j = 0; j < n_states; j++) {
			obs_prob += fwd[t][j] * bwd[t][j];
		}

		for (int i = 0; i < n_states; i++) {
			// joint probability of being in state `i` at time `t` and
			// observing the sequence `observation`
			double joint_prob = fwd[t][i] * bwd[t][i];

			// straightforward derivation from Bayes' theorem
			gamma[i][t] = joint_prob / obs_prob;
		}
	}

}


/* update the transition and emission matrices based on the expectations
 * calculated on the expectation step */
void maximization_step(Matrix_v& transition, Matrix_v& emission, vector<double>& init_prob,
		const vector<int>& observation, const Matrix_v& gamma, const vector<Matrix_v>& ksi) {

	int n_emissions = emission.size();
	int n_states = transition.size();
	int T = observation.size();

	for (int i = 0; i < n_states; i++) {
		/* expected number of transitions from state `i` */
		double n_trans_i = 0.0;

		for (int t = 0; t < T-1; t++) {
			for (int k = 0; k < n_states; k++) {
				n_trans_i += ksi[t][i][k];
			}
		}

		for (int j = 0; j < n_states; j++) {
			/* expected number of transitions from state `i` to state `j` */
			double n_trans_ij = 0.0;

			for (int t = 0; t < T-1; t++) {
				n_trans_ij += ksi[t][i][j];
			}

			/* update transition matrix based on likelihood maximization */
			transition[i][j] = n_trans_ij / n_trans_i;
		}
	}


	for (int i = 0; i < n_states; i++) {
		// expected number of times of being in state `i`
		double n_obs_i = 0.0;
		for (int t = 0; t < T; t++) {
			n_obs_i += gamma[i][t];
		}

		for (int m = 0; m < n_emissions; m++) {
			// expected number of times being in state `i` and observing symbol `m`
			double n_obs_im = 0.0;
			for (int t = 0; t < T; t++) {
				if (observation[t] == m) {
					n_obs_im += gamma[i][t];
				}
			}

			// update emission matrix based on likelihood maximization
			emission[m][i] = n_obs_im / n_obs_i;
		}
	}


	for (int i = 0; i < n_states; i++) {
		init_prob[i] = gamma[i][0];
	}

}


Matrix_v forward(const Matrix_v& transition, const Matrix_v& emission,
	 const vector<double>& init_prob, const vector<int>& observation,
	 vector<double>& scale_c) {

	int n_states = transition.size();
	int T = observation.size();

	Matrix_v fwd(T, vector<double>(n_states));


	/* calculate initial probability */
	for (int i = 0; i < n_states; i++) {
		fwd[0][i] = init_prob[i] * emission[observation[0]][i];
		scale_c[0] += fwd[0][i];
	}

	scale_c[0] = 1.0 / scale_c[0];

	/* scaling scheme initialization step */
	for (int i = 0; i < n_states; i++) {
		fwd[0][i] *= scale_c[0];
	}


	/* employ dynamic programming to drop exponential complexity to O(T*(N^2)) */
	for (int t = 1; t < T; t++) {

		for (int j = 0; j < n_states; j++) {
			for (int k = 0; k < n_states; k++) {
				fwd[t][j] += fwd[t-1][k] * transition[k][j] * emission[observation[t]][j];
			}
			scale_c[t] += fwd[t][j];
		}
		scale_c[t] = (1.0 / scale_c[t]);

		for (int j = 0; j < n_states; j++) {
			fwd[t][j] *= scale_c[t];
		}
	}

	return fwd;
}


/* operations = {3 * (T-1) * (N^2)} + N */
Matrix_v backward(const Matrix_v& transition, const Matrix_v& emission,
		const vector<double>& init_prob, const vector<int>& observation,
		const vector<double>& scale_c) {

	int n_states = transition.size();
	int T = observation.size();

	Matrix_v bwd(T, vector<double>(n_states));

	/* initialization step (n_ops = N) */
	for (int i = 0; i < n_states; i++) {
		bwd[T-1][i] = scale_c[T-1];
	}

	/* recursion - DP step { n_ops = 3 * (T-1) * (N^2)) } */
	for (int t = T-2; t >= 0; t--) {
		for (int i = 0; i < n_states; i++) {
			for (int j = 0; j < n_states; j++) {
				bwd[t][i] += transition[i][j] * emission[observation[t+1]][j] * bwd[t+1][j];
			}
			bwd[t][i] *= scale_c[t];
		}
	}

	return bwd;
}
