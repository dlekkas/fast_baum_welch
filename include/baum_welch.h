#ifndef BAUM_WELCH_H
#define BAUM_WELCH_H


#include <vector>
#include "hmm.h"

using IVector = std::vector<int>;

class BaumWelch {

	public:

		virtual void Load(HMM& hmm, IVector obs_seq) {};

		virtual void operator()() = 0; //{ std::cout << "Swsth" << std::endl;};

		virtual HMM GetHMM() { return HMM(0,0); };

		virtual ~BaumWelch() {};

};


#endif
