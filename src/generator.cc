#include <random>
#include <iostream>
#include <cassert>

#include "../include/generator.h"

using namespace std;


vector<double> dirichlet_sample(int n_samples, const vector<double>& alpha, double beta) {

	assert(n_samples <= (int) alpha.size());

	std::random_device rd;
	std::mt19937 gen(rd());

	vector<std::gamma_distribution<>> dist(n_samples);
	for (auto i = 0; i < n_samples; i++) {
		dist[i] = std::gamma_distribution<>(alpha[i], beta);
	}

	vector<double> samples(n_samples);
	for (auto i = 0; i < n_samples; i++) {
		samples[i] = dist[i](gen);
	}

	double norm = std::accumulate(samples.begin(), samples.end(), 0.0);
	for (auto& x: samples) {
		x /= norm;
	}

	return samples;
}



vector<double> symmetric_dirichlet_sample(int n_samples, double concentrate, double beta) {

	vector<double> alpha(n_samples, concentrate);
	return dirichlet_sample(n_samples, alpha, beta);
}


vector<double> uniform_dist_sample(int n_samples) {
	return symmetric_dirichlet_sample(n_samples, 1.0, 1.0);
}


vector<int> weighted_emission_sample(int n_samples, const vector<int>& weights) {

	std::random_device rd;
	std::mt19937 gen(rd());
	std::discrete_distribution<> dist(weights.begin(), weights.end());

	vector<int> res(n_samples);
	for (auto i = 0; i < n_samples; i++) {
		res[i] = dist(gen);
	}

	return res;
}

vector<int> uniform_emission_sample(int n_samples) {
	vector<double> weights(n_samples, 1.0);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::discrete_distribution<> dist(weights.begin(), weights.end());

	vector<int> res(n_samples);
	for (auto i = 0; i < n_samples; i++) {
		res[i] = dist(gen);
	}

	return res;
}
