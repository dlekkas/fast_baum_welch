//#include <mkl.h>
#include <mm_malloc.h>
#include "../include/tsc_x86.h"
#include "../include/benchmark.h"
#include <iostream>

#define data_size 1024*1024*128
#define n_iter 10
#define n_runs 10
#define alignment 64

using namespace std;

int main () {

    double* x = (double*)_mm_malloc(data_size*sizeof(double), alignment);
    double* y = (double*)_mm_malloc(data_size*sizeof(double), alignment);

    std::cout << "Start measurements" << std::endl;
    std::vector<double> cycles_list;

	init_tsc();

    for (int i=0; i<n_iter; i++) {

        uint64_t start = start_tsc();
		for (int j = 0; j < n_runs; j++) {
			//cblas_dcopy(size, x, 1, y, 1);
            for (int k=0; k<data_size; k++)
                x[k]=y[k];
		}
		uint64_t end = stop_tsc();
		uint64_t cycles = (end - start) / (double) n_runs;
        cycles_list.emplace_back(cycles);

    }

    Benchmark bench(cycles_list, "Memory bandwidth", "cycles", 0, 0, 0, 0);
    bench.BeautyPrint(cout);

    double b = data_size*8/ bench.stats.median;
    cout << "Measured Bandwidth: " << b << "bytes/cycle" << endl;



}

