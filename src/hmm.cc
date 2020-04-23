#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>

#include "../include/hmm.h"


void HMM::initialize_vectors(const std::string& input_file) {

    std::ifstream ifs {input_file};
    std::string line;

    ifs >> N >> M;

	std::cout << "N: " << N << std::endl;
	std::cout << "M: " << M << std::endl;

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
		std::getline(ifs, line);
		std::istringstream buf(line);
		std::copy(std::istream_iterator<double>(buf),
				  std::istream_iterator<double>(), A[i]);
	}

	std::getline(ifs, line);

    B = new double*[M];
	for (auto i = 0; i < M; i++) {
        B[i] = new double[N];
		std::getline(ifs, line);
		std::istringstream buf(line);
		std::copy(std::istream_iterator<double>(buf),
				  std::istream_iterator<double>(), B[i]);
	}

}


HMM::~HMM() {
	for (int i = 0; i < M; i++) {
		delete[] A[i];
	}
	for (int i = 0; i < N; i++) {
		delete[] B[i];
	}

	delete[] A;
	delete[] B;
}
