#include <iostream>
#include <vector>

#define MAX_ITER 10

using namespace std;

using Matrix_v = vector<vector<double>>;

Matrix_v forward(const Matrix_v&, const Matrix_v&, const vector<int>&,
		const vector<double>&);

Matrix_v backward(const Matrix_v&,  const Matrix_v&, const vector<int>&,
		const vector<double>&);

void maximization_step(Matrix_v&, Matrix_v&, const vector<int>&,
		const Matrix_v&, const vector<Matrix_v>&);

void expectation_step(const Matrix_v&, const Matrix_v&, const Matrix_v&,
		const Matrix_v&, const vector<int>&, Matrix_v&, vector<Matrix_v>&);




void baum_welch(Matrix_v& transition, Matrix_v& emission, vector<double>& init_prob,
		const vector<int>& observation) {

	int T = observation.size();
	int n_states = transition.size();

	/* probability of being in state `i` at time `t` and state `j` at time `t+1`,
	 * given the observation sequence and HMM model */
	vector<Matrix_v> ksi(T, Matrix_v(n_states, vector<double>(n_states)));

	/* probability of being in state `i` at time `t` */
	Matrix_v gamma(n_states, vector<double>(T));

	for (int i = 0; i < MAX_ITER; i++) {
		Matrix_v fwd = forward(transition, emission, observation, init_prob);
		Matrix_v bwd = backward(transition, emission, observation, init_prob);

		/* EXPECTATION step */
		expectation_step(transition, emission, fwd, bwd, observation, gamma, ksi);

		/* MAXIMIZATION step */
		maximization_step(transition, emission, observation, gamma, ksi);
	}

}


/* calculate the expected state occupancy count and the expected
 * state transition count based on the HMM parameters of previous step */
void expectation_step(const Matrix_v& transition, const Matrix_v& emission,
		const Matrix_v& fwd, const Matrix_v& bwd, const vector<int>& observation,
		Matrix_v& gamma, vector<Matrix_v>& ksi) {

	int n_states = transition.size();
	int T = observation.size();

	for (int t = 0; t < T; t++) {
		/* probability of observing sequence `observation` given the HMM model */
		double obs_prob = 0.0;
		for (int j = 0; j < n_states; j++) {
			obs_prob += fwd[t][j] * bwd[t][j];
		}

		for (int i = 0; i < n_states; i++) {
			for (int j = 0; j < n_states; j++) {
				/* probability of being in state `i` at time `t`, at state `j`
				 * at time `t+1` and observing sequence `observation` */
				double joint_prob = fwd[t][i] * transition[i][j] *
					emission[j][observation[t+1]] * bwd[t+1][j];

				/* straightforward derivation from Bayes' theorem */
				ksi[t][i][j] = joint_prob / obs_prob;
			}
		}
	}


	for (int i = 0; i < n_states; i++) {
		/* probability of observing sequence `observation` given the HMM model */
		double obs_prob = 0.0;
		for (int t = 0; t < T; t++) {
			obs_prob += fwd[t][i] * bwd[t][i];
		}

		for (int t = 0; t < T; t++) {
			/* joint probability of being in state `i` at time `t` and
			 * observing the sequence `observation` */
			double joint_prob = fwd[t][i] * bwd[t][i];

			/* straightforward derivation from Bayes' theorem */
			gamma[i][t] = joint_prob / obs_prob;
		}
	}
}


/* update the transition and emission matrices based on the expectations
 * calculated on the expectation step */
void maximization_step(Matrix_v& transition, Matrix_v& emission, const vector<int>& observation,
		const Matrix_v& gamma, const vector<Matrix_v>& ksi) {

	int n_emissions = emission[0].size();
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
		/* expected number of times of being in state `i` */
		double n_obs_i = 0.0;
		for (int t = 0; t < T; t++) {
			n_obs_i += gamma[i][t];
		}

		for (int m = 0; m < n_emissions; m++) {
			/* expected number of times being in state `i` and observing symbol `m` */
			double n_obs_im = 0.0;
			for (int t = 0; t < T; t++) {
				if (observation[t] == m) {
					n_obs_im += gamma[i][t];
				}
			}

			/* update emission matrix based on likelihood maximization */
			emission[i][m] = n_obs_im / n_obs_i;
		}
	}
}


/* operations = (3 * T * (N^2)) + N */
Matrix_v forward(const Matrix_v& transition, const Matrix_v& emission,
	 const vector<double>& init_prob, const vector<int>& observation) {

	int n_states = transition.size();
	int T = observation.size();

	Matrix_v fwd(T, vector<double>(n_states));

	/* calculate initial probability */
	/* operations = N */
	for (int i = 0; i < n_states; i++) {
		fwd[0][i] = init_prob[i] * emission[i][observation[0]];
	}

	/* employ dynamic programming to drop exponential complexity to O(T*(N^2)) */
	/* operations = 3 * T * (N^2) */
	for (int t = 1; t < T; t++) {
		for (int j = 0; j < n_states; j++) {
			for (int k = 0; k < n_states; k++) {
				fwd[t][j] += fwd[t-1][k] * transition[k][j] * emission[j][observation[t]];
			}
		}
	}

	return fwd;
}


/* operations = {3 * (T-1) * (N^2)} + N */
Matrix_v backward(const Matrix_v& transition, const Matrix_v& emission,
		const vector<double>& init_prob, const vector<int>& observation) {

	int n_states = transition.size();
	int T = observation.size();

	Matrix_v bwd(T, vector<double>(n_states));

	/* initialization step (n_ops = N) */
	for (int i = 0; i < n_states; i++) {
		bwd[i][T-1] = 1.0;
	}

	/* recursion - DP step { n_ops = 3 * (T-1) * (N^2)) } */
	for (int i = 0; i < n_states; i++) {
		for (int t = 0; t < T-1; t++) {
			for (int j = 0; j < n_states; j++) {
				bwd[i][t] += transition[i][j] * emission[j][observation[t+1]] * bwd[j][t+1];
			}
		}
	}

	return bwd;
}



