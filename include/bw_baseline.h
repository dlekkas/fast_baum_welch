#ifndef BW_BASELINE_H
#define BW_BASELINE_H

#include <vector>
#include "hmm.h"


using IVector = std::vector<int>;

#define MAX_ITERATIONS 5


class BaumWelch {

	public:

		virtual void Load(HMM& hmm, IVector obs_seq) {};

		virtual void operator()() = 0; //{ std::cout << "Swsth" << std::endl;};

		virtual HMM GetHMM() { return HMM(0,0); };

		virtual ~BaumWelch() {};

};

class BaumWelchCpp: public BaumWelch {

	public:

		BaumWelchCpp() {};

		void Load(HMM& hmm, IVector obs_seq) {
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

		DVector init_prob;

		IVector observation;

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
