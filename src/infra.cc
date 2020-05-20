#include <chrono>
#include <iostream>
#include <fstream>
#include <numeric>
#include <iterator>
#include <algorithm>

#include "../include/tsc_x86.h"
#include "../include/benchmark.h"
#include "../include/infra.h"
#include "../include/baum_welch.h"


// #define MEASURE_FLOPS

#ifdef MEASURE_FLOPS
#include "../../pcm/cpucounters.h"
#include "../../pcm/cpuasynchcounter.h"
#include "../../pcm/topology.h"
#endif

void perf_test_rdtscp(const std::string& impl_tag, BaumWelch* impl, int M, int N, int S,
		int n_runs, int n_iter, std::ostream& xout, bool to_CSV, std::string out_file) {

	std::vector<int> observations = uniform_emission_sample(S, N);
	std::vector<double> cycles_list;

	init_tsc();

	HMM base_model(M, N);
	for (auto i = 0; i < n_iter; i++) {
		impl->Load(base_model, observations);

		uint64_t start = start_tsc();
		for (auto j = 0; j < n_runs; j++) {
			(*impl)();
		}
		uint64_t end = stop_tsc();
		uint64_t cycles = (end - start) / (double) n_runs;

		cycles_list.emplace_back(cycles);
	}

	Benchmark bench(cycles_list, impl_tag, "cycles", N, M, observations.size(), MAX_ITERATIONS);

	if (to_CSV) {
		bench.CSVPrint(out_file);
	} else {
		bench.CompactPrint(xout);
	}

}


void perf_test_chrono(const std::string& impl_tag, BaumWelch* impl, int M, int N, int S,
		int n_runs, int n_iter, std::ostream& xout, bool to_CSV, std::string out_file) {

	std::vector<int> observations = uniform_emission_sample(S, N);
	std::vector<double> time_list;

	HMM base_model(M, N);
	for (auto i = 0; i < n_iter; i++) {
		impl->Load(base_model, observations);

		auto begin = std::chrono::steady_clock::now();
		for (auto j = 0; j < n_runs; j++) {
			(*impl)();
		}
		auto end = std::chrono::steady_clock::now();
		auto duration_us = std::chrono::duration_cast
			<std::chrono::milliseconds>(end - begin).count() / n_runs;

		time_list.emplace_back(duration_us);
	}

	Benchmark bench(time_list, impl_tag, "msec", base_model.N, base_model.M,
			observations.size(), MAX_ITERATIONS);

	if (to_CSV) {
		bench.CSVPrint(out_file);
	} else {
		bench.CompactPrint(xout);
	}
}

#ifdef MEASURE_FLOPS
void perf_test_rdtscp_and_flops(const std::string& impl_tag, BaumWelch* impl, int M, int N, int S,
		int n_runs, int n_iter, std::ostream& xout, bool to_CSV, std::string out_file) {

	std::vector<int> observations = uniform_emission_sample(S, N);
	std::vector<double> cycles_list;
	std::vector<uint64_t> flops_list;

	init_tsc();

	int NB_EVENTS = 3;
	PCM::CustomCoreEventDescription events[NB_EVENTS];
	double values[NB_EVENTS];
	uint64_t multiplier[NB_EVENTS];

	// cpu/umask=0x01,event=0xC7,name=FP_ARITH_INST_RETIRED.SCALAR_DOUBLE/       FLOPS x1
	events[0].umask_value  = 0x01;
	events[0].event_number = 0xc7;
	multiplier[0] = 1;
	// cpu/umask=0x04,event=0xC7,name=FP_ARITH_INST_RETIRED.128B_PACKED_DOUBLE/  FLOPS x2
	events[1].umask_value  = 0x04;
	events[1].event_number = 0xc7;
	multiplier[1] = 2;
	// cpu/umask=0x10,event=0xC7,name=FP_ARITH_INST_RETIRED.256B_PACKED_DOUBLE/  FLOPS x4
	events[2].umask_value  = 0x10;
	events[2].event_number = 0xc7;
	multiplier[2] = 4;

	HMM base_model(M, N);
	for (auto i = 0; i < n_iter; i++) {
		impl->Load(base_model, observations);

		PCM * m = PCM::getInstance();
		m->disableJKTWorkaround();
		m->resetPMU();
		if (m->program(PCM::CUSTOM_CORE_EVENTS,&events) != PCM::Success)
			std::cout << "WARNING!!!!: Something went wrong with PCM." << std::endl;

		SystemCounterState pcm_before = getSystemCounterState();
		uint64_t start = start_tsc();
		for (auto j = 0; j < n_runs; j++) {
			(*impl)();
		}
		uint64_t end = stop_tsc();
		SystemCounterState pcm_after = getSystemCounterState();
		uint64_t cycles = (end - start) / (double) n_runs;

		uint64_t flops = 0;
		for (int i = 0; i < NB_EVENTS; i++) {
			uint64 value = getNumberOfCustomEvents(i, pcm_before, pcm_after);
			// printf("Event %0d: 0x%02x 0x%02x: %lld\n", i+1, events[i].event_number, events[i].umask_value, value);
			flops += value * multiplier[i];
		}

		cycles_list.emplace_back(cycles);
		flops_list.emplace_back(flops);
	}

	Benchmark bench(cycles_list, flops_list, impl_tag, "cycles", N, M, observations.size(), MAX_ITERATIONS);

	if (to_CSV) {
		bench.CSVPrint(out_file);
	} else {
		bench.CompactPrint(xout);
	}

}
#endif


bool IsValidImpl(BaumWelch* impl) {
	BaumWelchCppBaseline base_impl;
	int n = 64, m = 64, o = 128;

	std::vector<int> obs = uniform_emission_sample(o, n);
	HMM base_model(m, n);

	impl->Load(base_model, obs);
	(*impl)();
	HMM test_hmm = impl->GetHMM();

	base_impl.Load(base_model, obs);
	base_impl();
	HMM base_hmm = base_impl.GetHMM();

	return base_hmm.IsSimilar(test_hmm);
}
