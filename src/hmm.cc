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

    int i;

    std::string tmp;
	std::getline(ifs, tmp);
    std::getline(ifs, tmp);
    std::getline(ifs, tmp);

    std::istringstream buf(tmp);
    std::vector<double> line { std::istream_iterator<double>(buf), std::istream_iterator<double>()};
	for (vector<double>::iterator it = line.begin(); it != line.end(); it++) {
		pi.push_back(*it);
	}
    std::getline(ifs, tmp);

    for (i=0; i<M; i++) {
        std::getline(ifs, tmp);
        std::istringstream buf(tmp);
        std::vector<double> line { std::istream_iterator<double>(buf), std::istream_iterator<double>()};
		std::vector<double> v;
		for (vector<double>::iterator it = line.begin(); it != line.end(); it++) {
		    v.push_back(*it);
	    }
        A.push_back(v);
    }

    std::getline(ifs, tmp);

    for (i=0; i<M; i++) {
        std::getline(ifs, tmp);
        std::istringstream buf(tmp);
        std::vector<double> line { std::istream_iterator<double>(buf), std::istream_iterator<double>()};
		std::vector<double> v;
		for (vector<double>::iterator it = line.begin(); it != line.end(); it++) {
		    v.push_back(*it);
	    }
        B.push_back(v);
    }
}
