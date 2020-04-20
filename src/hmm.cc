#include <vector>
#include <iostream>
#include <random>
#include <fstream> 
#include <sstream>
#include <iterator>

#include "../include/hmm.h"

using namespace std;

void HMM::initialize_vectors(char* input_file) {

    std::ifstream ifs;
	ifs.open(input_file, std::ios_base::in);
	if (ifs.fail()) {
		std::cerr << "Error opening file " << input_file << std::endl;
		std::exit(2);
	}

    ifs >> N >> M;

    cout << "N: " << N << endl;
    cout << "M: " << M << endl;

    int i=0, j=0;

    std::string tmp;
	std::getline(ifs, tmp);
    std::getline(ifs, tmp);
    std::getline(ifs, tmp);

    pi = (double*)malloc(M*sizeof(double));
    std::istringstream buf(tmp);
    std::vector<double> line { std::istream_iterator<double>(buf), std::istream_iterator<double>()};
	for (vector<double>::iterator it = line.begin(); it != line.end(); it++) {
		pi[i]=*it;
        i++;
	}
    std::getline(ifs, tmp);

    A = (double**)malloc(M*sizeof(double*));
    for (i=0; i<M; i++) {
        A[i] = (double*)malloc(M*sizeof(double));
        j=0;
        std::getline(ifs, tmp);
        std::istringstream buf(tmp);
        std::vector<double> line { std::istream_iterator<double>(buf), std::istream_iterator<double>()};
		for (vector<double>::iterator it = line.begin(); it != line.end(); it++) {
            A[i][j] = *it;
            j++;
	    }
    }

    std::getline(ifs, tmp);

    B = (double**)malloc(M*sizeof(double*));
    for (i=0; i<M; i++) {
        j=0;
        B[i] = (double*)malloc(N*sizeof(double));
        std::getline(ifs, tmp);
        std::istringstream buf(tmp);
        std::vector<double> line { std::istream_iterator<double>(buf), std::istream_iterator<double>()};
		for (vector<double>::iterator it = line.begin(); it != line.end(); it++) {
		    B[i][j] = *it;
            j++;
	    }
    }
}
