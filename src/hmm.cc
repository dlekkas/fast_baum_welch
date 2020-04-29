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
		emission.emplace_back(symmetric_dirichlet_sample(N));
	}

	alloc_mem();
	for (size_t i = 0; i < transition.size(); i++) {
		std::copy(transition[i].begin(), transition[i].end(), A[i]);
		std::copy(emission[i].begin(), emission[i].end(), B[i]);
	}
}


void HMM::InitParamsCustom(const Matrix_v& trans, const Matrix_v& emis,
		const std::vector<double>& init_prob) {
	transition = trans;
	emission = emis;
	pi = init_prob;

	alloc_mem();
	for (size_t i = 0; i < trans.size(); i++) {
		std::copy(transition[i].begin(), transition[i].end(), A[i]);
		std::copy(emission[i].begin(), emission[i].end(), B[i]);
	}
}


void HMM::alloc_mem() {
	for (auto i = 0; i < M; i++) {
		A[i] = new double[M];
		B[i] = new double[N];
	}
}


HMM::HMM(int states, int emissions):
	M(states),
	N(emissions),
	A(new double*[states]),
	B(new double*[states])
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
	B(new double*[hmm.M]),
	transition(hmm.transition),
	emission(hmm.emission)
{
	alloc_mem();
	std::copy(hmm.A, hmm.A + (M*M), A);
	std::copy(hmm.B, hmm.B + (M*N), B);
}



HMM::~HMM() {
	if (A != nullptr) {
		for (int i = 0; i < M; i++) {
			delete[] A[i];
		}
		delete[] A;
	}

	if (B != nullptr) {
		for (int i = 0; i < M; i++) {
			delete[] B[i];
		}
		delete[] B;
	}
}
