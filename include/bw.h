#ifndef BW_H
#define BW_H

#include <vector>
#include <unordered_map>
#include <iostream>
#include <random>

#include "hmm.h"

#define TH 1e-4
#define MAX_ITERATIONS 5

class BW {

	public:

        HMM *hmm;
		
		BW(HMM *init_model, std::vector<int> seq, int t): hmm(init_model), threshold(TH), max_iterations(MAX_ITERATIONS), observation_seq(seq), T(t) {};

		~BW() {};

        void forward_backward(double** forward, double** backward);

        bool update_and_check(double** forward, double** backward);

        void run_bw();

        // convergence threshold
        double threshold;

        // maximum number of BW iterations
        int max_iterations;

        // observation sequence -- this should be a set of observation sequences
        std::vector<int> observation_seq;
        
        // observation sequence length
        int T;

};


#endif 			
