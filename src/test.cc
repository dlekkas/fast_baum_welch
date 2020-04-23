#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>

#include "../include/bw.h"


int main(int argc, char** argv) {

    if (argc < 3) {
		std::cout << "Usage: ./baum_welch <initialization file> \
			<observations file>" << std::endl;
		std::exit(1);
    }

    HMM model(argv[1]);

    std::vector<int> observations;
    std::ifstream ifs(argv[2]);
	std::copy(std::istream_iterator<double>(ifs), std::istream_iterator<double>(),
			  std::back_inserter(observations));

    run_bw(model.M, model.N, observations.size(), observations.data(), model.pi.data(), model.A, model.B);

}
