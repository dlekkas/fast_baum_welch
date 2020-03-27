#ifndef BW_H
#define BW_H

#include <vector>
#include <unordered_map>
#include <iostream>
#include <random>

#include "hmm.h"

#define TH 0.001

class BW {

	public:

        HMM *hmm;
		
		BW(HMM *init_model, int* seq): hmm(init_model), threshold(TH), observation_seq(seq) {};

		~BW() {};

        void forward_backward(double** a, double** b);

        void update(double** a, double** b);

        void check_convergence();

	private:

        // convergence threshold
        double threshold;

        // observation sequence -- this should be a set of observation sequences
        int* observation_seq;

};


#endif 			
