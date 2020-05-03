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
};


class Benchmark {

	public:

		int N; int M; int O;

		std::vector<double> measurements;
		std::string impl_tag;
		std::string metric_tag;

		Statistics stats;

		Benchmark(const std::vector<double>& values, const std::string& i_tag,
				const std::string& m_tag, int n, int m, int o);

		Benchmark(const Benchmark& bench) = default;

		~Benchmark() = default;

		Statistics CalculateStatistics();

		void BeautyPrint(std::ostream& os);

		void CompactPrint(std::ostream& os);

		void CSVPrint(const std::string& file);


};




#endif
