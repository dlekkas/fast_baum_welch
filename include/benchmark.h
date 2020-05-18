#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <vector>
#include <string>

struct Statistics {
	double mean;
	double median;
	double min_v;
	double max_v;
	double variance;
	double confidence_interval_high;
	double confidence_interval_low;
};


class Benchmark {

	public:

		int N; int M; int O;
		int bw_iterations;

		std::vector<double> measurements;
		std::string impl_tag;
		std::string metric_tag;

		Statistics stats;

		Benchmark(const std::vector<double>& values, const std::string& i_tag,
				const std::string& m_tag, int n, int m, int o, int bw_iters);

		Benchmark(const Benchmark& bench) = default;

		~Benchmark() = default;

		Statistics CalculateStatistics();

		void BeautyPrint(std::ostream& os);

		void CompactPrint(std::ostream& os);

		void CSVPrint(const std::string& file);

};

#endif
