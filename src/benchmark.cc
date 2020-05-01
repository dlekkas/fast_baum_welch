#include <algorithm>
#include <numeric>
#include <iostream>

#include "benchmark.h"


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


void Benchmark::BeautyPrint(std::ostream& os) {

	os << "--------------- " << impl_tag << "--------------- " << std::endl;
	os << "[" << metric_tag << "] (MEAN): "<< stats.mean << std::endl;
	os << "[" << metric_tag << "] (MEDIAN): "<< stats.median << std::endl;
	os << "[" << metric_tag << "] (VARIANCE): "<< stats.median << std::endl;
	os << "[" << metric_tag << "] (MIN): " << stats.min_v << std::endl;
	os << "[" << metric_tag << "] (MIN): " << stats.max_v << std::endl;

}


Benchmark::Benchmark(const std::vector<double>& values, const std::string& i_tag,
		const std::string& m_tag, int n, int m, int o):
	N(n), M(m), O(o),
	measurements(values),
	impl_tag(i_tag),
	metric_tag(m_tag)
{
	CalculateStatistics();
}
