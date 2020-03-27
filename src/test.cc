#include <vector>
#include <iostream>
#include <random>
#include <fstream> 
#include <sstream>
#include <iterator>

#include "../include/bw.h"

using namespace std;

int main(int argc, char** argv) {

    if (argc < 3) {
        cout << "Usage: ./baum_welch <initialization file> <observations file>" << endl;
        return 0;
    }    

    HMM* model = new HMM();
    model->initialize_vectors(argv[1]);
    
    std::ifstream ifs;
	ifs.open(argv[2], std::ios_base::in);
	if (ifs.fail()) {
		std::cerr << "Error opening file " << argv[2] << std::endl;
		std::exit(2);
	}

    std::vector<int> v;
    std::string tmp;
	std::getline(ifs, tmp);
    std::istringstream buf(tmp);
    std::vector<double> line { std::istream_iterator<double>(buf), std::istream_iterator<double>()};
	for (vector<double>::iterator it = line.begin(); it != line.end(); it++) {
		v.push_back(*it);
	}

    BW* baum_welch = new BW(model, v, v.size());
    baum_welch->run_bw();

    return 1;
}
