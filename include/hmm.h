#ifndef HMM_H
#define HMM_H

#include <string>
#include <vector>
#include <iostream>

#include "generator.h"

using Matrix = double**;

using Matrix_v = std::vector<std::vector<double>>;
using DVector = std::vector<double>;


class HMM {

	public:

		HMM(int states, int emissions);

		HMM(const HMM& hmm);

		HMM(Matrix_v trans, Matrix_v emis, DVector init_prob);


		HMM(const std::string& input_file);

		~HMM();


		void InitParamsCustom(const Matrix_v& trans, const Matrix_v& emis,
				const std::vector<double>& init_prob);

		//void InitParamsFromFile(const std::string& input_file);

		void InitParamsRandom();

		bool IsSimilar(const HMM& hmm, const double eps = 1e-4);


        // number of states
        int M;

        // number of observations
		int N;

		Matrix_v transition;

		Matrix_v emission;

        // initial state vector
		DVector pi;

};


#endif
