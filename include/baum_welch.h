#ifndef BAUM_WELCH_H
#define BAUM_WELCH_H

#include <vector>
#include "hmm.h"



#define MAX_ITERATIONS 5


class BaumWelch {

	public:

		virtual void Load(HMM& hmm, std::vector<int>& obs_seq) = 0;

		virtual void operator()() = 0;

		virtual HMM GetHMM() = 0;

		virtual ~BaumWelch() = default;

};


class BaumWelchC: public BaumWelch {

	public:

		BaumWelchC() = default;

		virtual void operator()() = 0;

		virtual void Load(HMM& hmm, std::vector<int>& obs_seq);

		HMM GetHMM();

		virtual ~BaumWelchC();


	protected:

		int M; int N; int T;

		double** A;

		double** B;

		double* pi;

		double** fwd;

		double** bwd;

		int* obs;

};

class BaumWelchCLoopUnroll: public BaumWelchC {
	public:
		void operator()();
};

class BaumWelchCLoopUnroll0: public BaumWelchC {
	public:
		void operator()();
};

class BaumWelchCVect: public BaumWelch {

	public:

		BaumWelchCVect() = default;

		virtual void operator()() = 0;

		virtual void Load(HMM& hmm, std::vector<int>& obs_seq);

		HMM GetHMM();

		virtual ~BaumWelchCVect();


	protected:

		int M; int N; int T;

		double* A;

		double* B;

		double* pi;

		double* fwd;

		double* bwd;

		int* obs;

};

class BaumWelchCVectBasic: public BaumWelchCVect {
	public:
		void operator()();
};

class BaumWelchCVectDim: public BaumWelchCVect {
	public:
		void operator()();
};





class BaumWelchCExtended: public BaumWelchC {

	public:
		void operator()() = 0;

		virtual void Load(HMM& hmm, std::vector<int>& obs_seq);

		virtual ~BaumWelchCExtended();

	protected:

		double** gamma;

		double*** chsi;

};


class BaumWelchCBasic: public BaumWelchCExtended {
	public:
		void operator()();
};

class BaumWelchCBasicOpts: public BaumWelchCExtended {
	public:
		void operator()();
};

class BaumWelchCOptsV2: public BaumWelchCExtended {
	public:
		void operator()();
};





class BaumWelchCpp: public BaumWelch {

	public:

		BaumWelchCpp() = default;

		void Load(HMM& hmm, std::vector<int>& obs_seq) {
			transition = hmm.transition;
			emission = hmm.emission;
			init_prob = hmm.pi;
			observation = obs_seq;
		}

		HMM GetHMM() {

			HMM res(transition, emission, init_prob);

			return res;
		}

		virtual void operator()() = 0;

		~BaumWelchCpp() {};


	protected:

		Matrix_v transition;

		Matrix_v emission;

		std::vector<double> init_prob;

		std::vector<int> observation;

};


class BaumWelchCppOpts: public BaumWelchCpp {
	public:
		void operator()();
};

class BaumWelchCppBaseline: public BaumWelchCpp {
	public:
		void operator()();
};




#endif
