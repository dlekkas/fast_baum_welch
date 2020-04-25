#ifndef HMM_H
#define HMM_H

#include <string>
#include <vector>

#include "generator.h"

using Matrix = double**;

using Matrix_v = std::vector<std::vector<double>>;

class HMM {

	public:

		HMM(): A(nullptr), B(nullptr) {};

		HMM(const std::string& input_file): A(nullptr), B(nullptr) {
			InitParamsFromFile(input_file);
		}

		~HMM();

		void InitParamsFromFile(const std::string& input_file);


		void InitParamsRandom();


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

		Matrix_v transition;

		Matrix_v emission;

};


#endif
