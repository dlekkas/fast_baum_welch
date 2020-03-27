#ifndef HMM_H
#define HMM_H


class HMM {

	public:
		
		HMM() {};

		~HMM() {};

		void initialize_vectors(char* input_file);

        // number of observations
		int N;

        // number of states
        int M;

        // initial state vector
        std::vector<double> pi;

        // state transition matrix
        std::vector<std::vector<double>> A;

        // emission matrix
        std::vector<std::vector<double>> B;


};


#endif 			
