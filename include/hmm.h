#ifndef HMM_H
#define HMM_H

#include <string>


class HMM {

	public:

		HMM() {};

		HMM(const std::string& input_file) {
			initialize_vectors(input_file);
		}

		~HMM() {};

		void initialize_vectors(const std::string& input_file);

        // number of observations
		int N;

        // number of states
        int M;

        // initial state vector
        double* pi;

        // state transition matrix
        double** A;

        // emission matrix
        double** B;

};


#endif
