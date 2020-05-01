#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>

#include "../include/bw.h"
#include "../include/infra.h"
#include "../include/hmm.h"
#include "../include/validation.h"
#include "../include/generator.h"

#define N_ITERATIONS 10
#define N_RUNS 1
#define SEQ_LEN 50


int main(int argc, char** argv) {

	/*
	 *  This is how we would get initialization parameters + observations from file
		std::string init_file { "init_params.txt" };
		std::string obs_file { "observations.txt" };
		perf_test_chrono(init_file, obs_file, &run_bw, N_RUNS, N_ITERATIONS, std::cout);
		perf_test_rdtscp(init_file, obs_file, &run_bw, N_RUNS, N_ITERATIONS, std::cout);
	*/

	int n_states = 100;
	int n_emissions = 100;

	//VALIDATION PART - MANOLIS

	std::vector<int> observations = uniform_emission_sample(SEQ_LEN, n_emissions);

	HMM model_old(n_states, n_emissions);
	HMM model_base(model_old);
	// C++ like implementation as basic
	baum_welch(model_base.transition, model_base.emission, model_base.pi, observations);

	for (size_t i = 0; i < model_base.transition.size(); i++) {
		std::copy(model_base.transition[i].begin(), model_base.transition[i].end(), model_base.A[i]);
		std::copy(model_base.emission[i].begin(), model_base.emission[i].end(), model_base.B[i]);
	}

    for (int j = 0; j < 1; j++) { // BE CAREFUL! Only one function for the moment, we compare the basic with itself!

		HMM model(model_old);
  	    //comp_func f = userFuncs[i];
		/*
		baum_welch(model.transition, model.emission, model.pi, observations);
		for (size_t i = 0; i < model.transition.size(); i++) {
			std::copy(model.transition[i].begin(), model.transition[i].end(), model.A[i]);
			std::copy(model.emission[i].begin(), model.emission[i].end(), model.B[i]);
		}
		*/
		run_bw(model.M, model.N, observations.size(), observations.data(), model.pi.data(), model.A, model.B);
		double error = compute_error(model_base.A, model_base.B, model_base.pi.data(), model.A, model.B, model.pi.data(), n_states, n_emissions);

        if (error > ERROR_BOUND) {
            std::cout << error << std::endl;
            //cout << "ERROR!!!!  the results for the " << i+1 << "th function are different to the previous" << std::endl;
			std::cout << "ERROR!!! The results of given function are different to basic implementation" << std::endl;
            return 1;
        }
    }
  // Add Garbage COllector

  //END OF VALIDATION PART

	// Measurements regarding the C++ baseline implementation
	perf_test_rdtscp("Baseline implementation", &baum_welch,
			n_states, n_emissions, N_RUNS, N_ITERATIONS, std::cout);
	perf_test_chrono("Baseline implementation", &baum_welch,
			n_states, n_emissions, N_RUNS, N_ITERATIONS, std::cout);

	// Measurements regarding the C-like implementation
	perf_test_rdtscp("Basic opts", &run_bw, n_states, n_emissions, N_RUNS, N_ITERATIONS, std::cout);
	perf_test_chrono("Basic opts", &run_bw, n_states, n_emissions, N_RUNS, N_ITERATIONS, std::cout);
}
