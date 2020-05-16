#ifndef INFRA_H
#define INFRA_H

#include <iostream>

#include "hmm.h"
#include "baum_welch.h"

void perf_test_rdtscp(const std::string& tag, BaumWelch* impl,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV = false,
		std::string out_file = "");

void perf_test_chrono(const std::string& tag, BaumWelch* impl,
		int M, int N, int S, int n_runs, int n_iter, std::ostream& xout, bool to_CSV = false,
		std::string out_file = "");

bool IsValidImpl(BaumWelch* impl);

#endif
