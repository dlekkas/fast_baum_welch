#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <math.h>

#include "../include/hmm.h"


/*
void HMM::InitParamsFromFile(const std::string& input_file) {

    std::ifstream ifs {input_file};
    std::string line;

    ifs >> N >> M;

	std::getline(ifs, line);
	std::getline(ifs, line);

    std::getline(ifs, line);
	std::istringstream buf {line};
	std::copy(std::istream_iterator<double>(buf),
			  std::istream_iterator<double>(), std::back_inserter(pi));

	std::getline(ifs, line);

    A = new double*[M];
	for (auto i = 0; i < M; i++) {
		A[i] = new double[M];
		std::getline(ifs, line); std::istringstream buf(line);
		std::copy(std::istream_iterator<double>(buf),
				  std::istream_iterator<double>(), A[i]);
	}

	std::getline(ifs, line);

	for (auto i = 0; i < M; i++) {
        B[i] = new double[N];
		std::getline(ifs, line); std::istringstream buf(line);
		std::copy(std::istream_iterator<double>(buf),
				  std::istream_iterator<double>(), B[i]);
	}

}
*/


void HMM::InitParamsRandom() {
	pi = symmetric_dirichlet_sample(M);

	transition.reserve(M);
	emission.reserve(M);
	for (auto i = 0; i < M; i++) {
		transition.emplace_back(symmetric_dirichlet_sample(M));
	}

	for (auto i = 0; i < N; i++) {
		emission.emplace_back(symmetric_dirichlet_sample(M));
	}
}


void HMM::InitParamsCustom(const Matrix_v& trans, const Matrix_v& emis,
		const std::vector<double>& init_prob) {
	transition = trans;
	emission = emis;
	pi = init_prob;
}




bool HMM::IsSimilar(const HMM& hmm, const double eps) {

	double max_err = 0.0;
	for (size_t i = 0; i < pi.size(); i++) {
		if (isnan(hmm.pi[i]) || isnan(pi[i]))
			return false;
		max_err = std::max(max_err, std::abs(pi[i] - hmm.pi[i]));
	}

	for (size_t i = 0; i < transition.size(); i++) {
		for (size_t j = 0; j < transition[0].size(); j++) {
			if (isnan(hmm.transition[i][j]) || isnan(transition[i][j]))
				return false;
			max_err = std::max(max_err, std::abs(transition[i][j] - hmm.transition[i][j]));
		}
	}

	for (size_t i = 0; i < emission.size(); i++) {
		for (size_t j = 0; j < emission[0].size(); j++) {
			if (isnan(hmm.emission[i][j]) || isnan(emission[i][j]))
				return false;
			max_err = std::max(max_err, std::abs(emission[i][j] - hmm.emission[i][j]));
		}
	}

	return max_err < eps;
}






HMM::HMM(int states, int emissions):
	M(states),
	N(emissions)
{
	InitParamsRandom();
}

HMM::HMM(const std::string& input_file):
	M(0), N(0) {
	//InitParamsFromFile(input_file);
}

HMM::HMM(Matrix_v trans, Matrix_v emis, DVector init_prob):
	transition(trans),
	emission(emis),
	pi(init_prob)
{}

HMM::HMM(const HMM& hmm):
	M(hmm.M),
	N(hmm.N),
	transition(hmm.transition),
	emission(hmm.emission),
	pi(hmm.pi)
{}

HMM::~HMM() {}
