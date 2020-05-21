
#include "../include/baum_welch.h"
#include "../include/hmm.h"
#include "../include/generator.h"



#define T 256
#define M 256
#define N 256

int main() {

	//ANNOTATE_DISABLE_COLLECTION_PUSH;
	HMM base_model(M, N);
	std::vector<int> observations = uniform_emission_sample(T, N);


	BaumWelchCppBaseline impl1;
	impl1.Load(base_model, observations);
	impl1();

	BaumWelchCppOpts impl2;
	impl2.Load(base_model, observations);
	impl2();

	BaumWelchCBasic impl3;
	impl3.Load(base_model, observations);
	impl3();

	BaumWelchCOptsV2 impl5;
	impl5.Load(base_model, observations);
	impl5();

	BaumWelchCOptsBlocking impl7;
	impl7.Load(base_model, observations);
	impl7();

	BaumWelchCLoopUnroll impl11;
	impl11.Load(base_model, observations);
	impl11();

	BaumWelchCVectBasic impl10;
	impl10.Load(base_model, observations);
	impl10();

	BaumWelchCVectDim2 impl12;
	impl12.Load(base_model, observations);
	impl12();

	BaumWelchCVectUnroll impl13;
	impl13.Load(base_model, observations);
	impl13();

	BaumWelchCVectOpt impl8;
	impl8.Load(base_model, observations);
	impl8();


	//ANNOTATE_DISABLE_COLLECTION_POP;


}




