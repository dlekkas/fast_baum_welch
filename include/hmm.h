#ifndef HMM_H
#define HMM_H

#include <string>
#include <vector>

using Matrix = double**;

class HMM {

	public:

		HMM() {};

		HMM(const std::string& input_file) {
			initialize_vectors(input_file);
		}

		~HMM();

		void initialize_vectors(const std::string& input_file);

        // number of observations
		int N;

        // number of states
        int M;

        // initial state vector
		std::vector<double> pi;

        // state transition matrix
		Matrix A;

        // emission matrix
        Matrix B;

};


#endif
