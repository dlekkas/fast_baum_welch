#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>

#include "benchmark.h"
#include "bw.h"

using namespace std;


Statistics Benchmark::CalculateStatistics() {
	int n = measurements.size();
	std::sort(measurements.begin(), measurements.end());

	stats.median = (n % 2 != 0) ? measurements[n/2] :
		(measurements[n/2] + measurements[n/2+1]) / 2;
	stats.mean = std::accumulate(measurements.begin(),
			measurements.end(), 0.0) / n;
	stats.min_v = measurements.front();
	stats.max_v = measurements.back();

	int sum = 0;
	for (auto x: measurements) {
		sum += (stats.mean-x)*(stats.mean-x);
	}
	stats.variance = sum / (n-1);

	return stats;
}


void Benchmark::BeautyPrint(ostream& os) {

	os << "--------------- " << impl_tag << "--------------- " << endl;
	os << "[" << metric_tag << "] (MEAN): "<< stats.mean << endl;
	os << "[" << metric_tag << "] (MEDIAN): "<< stats.median << endl;
	os << "[" << metric_tag << "] (VARIANCE): "<< stats.variance << endl;
	os << "[" << metric_tag << "] (MIN): " << stats.min_v << endl;
	os << "[" << metric_tag << "] (MAX): " << stats.max_v << endl;

}

void Benchmark::CSVPrint(const string& file) {
	ofstream ofs(file, ios_base::app);
	ofs << impl_tag << "," << metric_tag << "," << M << ","
		<< N << "," << O << "," << stats.mean << "," << bw_iterations << endl;
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
