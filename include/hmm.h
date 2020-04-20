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
        double* pi;

        // state transition matrix
        double** A;

        // emission matrix
        double** B;


};


#endif 			
