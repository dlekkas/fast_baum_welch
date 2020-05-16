#ifndef BW_BASELINE_H
#define BW_BASELINE_H

#include <vector>
#include "baum_welch.h"

using IVector = std::vector<int>;


class BaumWelchBaseline: public BaumWelch {

	public:

		BaumWelchBaseline() {};

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

		void operator()();

		~BaumWelchBaseline() {};


	private:

		Matrix_v transition;

		Matrix_v emission;

		DVector init_prob;

		IVector observation;

};

extern BaumWelchBaseline x;


#endif
