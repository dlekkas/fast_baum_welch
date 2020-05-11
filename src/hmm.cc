#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>

#include "../include/hmm.h"


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

    B = new double*[M];
	for (auto i = 0; i < M; i++) {
        B[i] = new double[N];
		std::getline(ifs, line); std::istringstream buf(line);
		std::copy(std::istream_iterator<double>(buf),
				  std::istream_iterator<double>(), B[i]);
	}

}


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

	alloc_mem();

	for (auto i = 0; i < N; i++) {
		std::copy(emission[i].begin(), emission[i].end(), B[i]);
	}


	for (auto i = 0; i < M; i++) {
		std::copy(transition[i].begin(), transition[i].end(), A[i]);
	}
}


void HMM::InitParamsCustom(const Matrix_v& trans, const Matrix_v& emis,
		const std::vector<double>& init_prob) {
	transition = trans;
	emission = emis;
	pi = init_prob;

	alloc_mem();
	for (auto i = 0; i < M; i++) {
		std::copy(transition[i].begin(), transition[i].end(), A[i]);
	}

	for (auto i = 0; i < N; i++) {
		std::copy(emission[i].begin(), emission[i].end(), B[i]);
	}
}



bool HMM::IsSimilar(const HMM& hmm, const double eps) {

	double max_err = 0.0;
	for (auto i = 0; i < M; i++) {
		max_err = std::max(max_err, std::abs(pi[i] - hmm.pi[i]));
	}

	for (auto i = 0; i < M; i++) {
		for (auto j = 0; j < M; j++) {
			max_err = std::max(max_err, std::abs(A[i][j] - hmm.A[i][j]));
		}
	}

	for (auto i = 0; i < N; i++) {
		for (auto j = 0; j < M; j++) {
			max_err = std::max(max_err, std::abs(B[i][j] - hmm.B[i][j]));
		}
	}

	return max_err < eps;
}




void HMM::alloc_mem() {
	for (auto i = 0; i < M; i++) {
		A[i] = new double[M];
	}

	for (auto i = 0; i < N; i++) {
		B[i] = new double[M];
	}
}


HMM::HMM(int states, int emissions):
	M(states),
	N(emissions),
	A(new double*[states]),
	B(new double*[emissions])
{
	InitParamsRandom();
}

HMM::HMM(const std::string& input_file):
	M(0), N(0), A(nullptr), B(nullptr) {
	InitParamsFromFile(input_file);
}


HMM::HMM(const HMM& hmm):
	M(hmm.M),
	N(hmm.N),
	pi(hmm.pi),
	A(new double*[hmm.M]),
	B(new double*[hmm.N]),
	transition(hmm.transition),
	emission(hmm.emission)
{
	alloc_mem();
  	for (auto i = 0; i < hmm.M; i++) {
		std::copy(transition[i].begin(), transition[i].end(), A[i]);
	}

  	for (auto i = 0; i < hmm.N; i++) {
		std::copy(emission[i].begin(), emission[i].end(), B[i]);
	}

}



HMM::~HMM() {
	if (A != nullptr) {
		for (int i = 0; i < M; i++) {
			delete[] A[i];
		}
		delete[] A;
	}



	if (B != nullptr) {
		for (int i = 0; i < N; i++) {
			delete[] B[i];
		}
		delete[] B;
	}
}
