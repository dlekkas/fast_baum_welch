#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>

#include "../include/bw.h"
#include "../include/benchmark.h"
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

	int n_states = 5;
	int n_emissions = 3;

	//VALIDATION PART - MANOLIS
  HMM model_old(n_states, n_emissions);
  HMM model_base(n_states, n_emissions);

	HMM model_1(n_states, n_emissions);
	model_1.InitParamsRandom();

	std::vector<int> observations = uniform_emission_sample(SEQ_LEN, n_emissions);

	model_old.InitParamsCustom(model_1.transition, model_1.emission, model_1.pi);
	// We use the C++ like implementation as Basic
	baum_welch(model_1.transition, model_1.emission, model_1.pi, observations);
	model_base.InitParamsCustom(model_1.transition, model_1.emission, model_1.pi);

  for (int j = 0; j < 1; j++) { // BE CAREFUL! Only one function for the moment, we compare the basic with itself!
		//std::vector<int> observations_2 = uniform_emission_sample(SEQ_LEN, n_emissions);
		model_1.InitParamsCustom(model_old.transition, model_old.emission, model_old.pi);
    //comp_func f = userFuncs[i];
  	baum_welch(model_1.transition, model_1.emission, model_1.pi, observations);

		for (size_t i = 0; i < model_1.transition.size(); i++) {
			std::copy(model_1.transition[i].begin(), model_1.transition[i].end(), model_1.A[i]);
			std::copy(model_1.emission[i].begin(), model_1.emission[i].end(), model_1.B[i]);
		}

    //double error = compute_error(model_base.transition, model_base.emission, model_base.pi, model_1.A, model_1.B, model_1.pi, M, N);
		double error = compute_error(model_base.A, model_base.B, model_base.pi.data(), model_1.A, model_1.B, model_1.pi.data(), n_states, n_emissions);

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
