#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <cmath>

#include "benchmark.h"
#include "bw.h"

using namespace std;


Statistics Benchmark::CalculateStatistics() {
	int n = measurements.size();
	std::sort(measurements.begin(), measurements.end());
	std::sort(flop_measurements.begin(), flop_measurements.end());

	stats.median = (n % 2 != 0) ? measurements[n/2] :
		(measurements[n/2] + measurements[n/2+1]) / 2;
	stats.mean = std::accumulate(measurements.begin(),
			measurements.end(), 0.0) / n;
	stats.min_v = measurements.front();
	stats.max_v = measurements.back();

	uint64_t sum = 0;
	for (auto x: measurements) {
		sum += (stats.mean-x)*(stats.mean-x);
	}
	stats.variance = sum / (n-1);

	// Calculation of 95% confidence interval according to
	// https://spcl.inf.ethz.ch/Teaching/2017-dphpc/hoefler-scientific-benchmarking.pdf; z(0.025) = 1.96
	int ci_low_rank  = floor((n - (1.96 * sqrt(n))) / 2);
	int ci_high_rank = ceil(1 + ((n + (1.96 * sqrt(n))) / 2));
	// check that they do not exceed the limits of array indexes
	ci_low_rank  = max(0, ci_low_rank - 1);
	ci_high_rank = min(n - 1, ci_high_rank - 1);

	stats.confidence_interval_low  = measurements[ci_low_rank];
	stats.confidence_interval_high = measurements[ci_high_rank];

	int n_flops = flop_measurements.size();
	if (n_flops != 0) {
		flop_stats.median = (n_flops % 2 != 0) ? (double)flop_measurements[n_flops/2] :
								((double)(flop_measurements[n_flops/2] + flop_measurements[n_flops/2+1])) / 2;
	}

#ifdef DEBUG_CI
	cout << "measurements:" << endl;
	for (auto &x : measurements) {
		cout << x << " ";
	}
	cout << "stats.confidence_interval_high: " << stats.confidence_interval_high << endl;
	cout << "stats.confidence_interval_low: " << stats.confidence_interval_low << endl;
	cout << "stats.median: " << stats.median << endl;

	cout << endl;
#endif

	return stats;
}


void Benchmark::BeautyPrint(ostream& os) {

	os << "--------------- " << impl_tag << "--------------- " << endl;
	os << "[" << metric_tag << "] (MEAN): "<< stats.mean << endl;
	os << "[" << metric_tag << "] (MEDIAN): "<< stats.median << endl;
	os << "[" << metric_tag << "] (VARIANCE): "<< stats.variance << endl;
	os << "[" << metric_tag << "] (MIN): " << stats.min_v << endl;
	os << "[" << metric_tag << "] (MAX): " << stats.max_v << endl;
	os << "[" << metric_tag << "] (CI_LOW_BOUND): " << stats.confidence_interval_low << endl;
	os << "[" << metric_tag << "] (CI_HIGH_BOUND): " << stats.confidence_interval_high << endl;

}

void Benchmark::CSVPrint(const string& file) {
	ofstream ofs(file, ios_base::app);
	ofs << impl_tag << "," << metric_tag << "," << M << ","
		<< N << "," << O << "," << stats.median << "," << bw_iterations << "," << stats.confidence_interval_low << ","
		<< stats.confidence_interval_high;

	if (flop_stats.median > 1) {
		ofs << "," << flop_stats.median;
	}
	ofs << endl;
}

void Benchmark::CompactPrint(ostream& os) {
	os << impl_tag << ": " << stats.median << " " << metric_tag << endl;
}


Benchmark::Benchmark(const vector<double>& values, const string& i_tag,
		const string& m_tag, int n, int m, int o, int bw_iters):
	N(n), M(m), O(o),
	bw_iterations(bw_iters),
	measurements(values),
	impl_tag(i_tag),
	metric_tag(m_tag)
{
	CalculateStatistics();
}

Benchmark::Benchmark(const vector<double>& values, const std::vector<uint64_t>& flop_values, const string& i_tag,
		const string& m_tag, int n, int m, int o, int bw_iters):
	N(n), M(m), O(o),
	bw_iterations(bw_iters),
	measurements(values),
	flop_measurements(flop_values),
	impl_tag(i_tag),
	metric_tag(m_tag)
{
	CalculateStatistics();
}
