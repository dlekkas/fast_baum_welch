#ifndef HMM_H
#define HMM_H

#include <string>
#include <vector>
#include <iostream>

#include "generator.h"

using Matrix = double**;

using Matrix_v = std::vector<std::vector<double>>;


class HMM {

	public:

		HMM(int states, int emissions);

		HMM(const HMM& hmm);

		HMM(const std::string& input_file);

		~HMM();


		void InitParamsCustom(const Matrix_v& trans, const Matrix_v& emis,
				const std::vector<double>& init_prob);

		void InitParamsFromFile(const std::string& input_file);

		void InitParamsRandom();

		bool IsSimilar(const HMM& hmm, const double eps = 1e-3);


        // number of states
        int M;

        // number of observations
		int N;

        // initial state vector
		std::vector<double> pi;

        // state transition matrix
		Matrix A;

        // emission matrix
        Matrix B;

		Matrix_v transition;

		Matrix_v emission;

	private:

		void alloc_mem();

};


#endif
