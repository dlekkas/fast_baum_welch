# Blazingly Fast Baum Welch for Hidden Markov Models

## Introduction
This repository contains the fastest publicly available Baum Welch implementation to our knowledge built in C/C++. Baum Welch is used to estimate the parameters of a Hidden Markov Model and has various practical applications ranging from speech recognition and cryptanalysis to algorithmic trading. The goal of this project was to explore a huge amount of optimizations in order to achieve the best attainable performance while targeting the *Skylake Micro-architecture*. This repo also contains the necessary benchmarking and validation infrastructure to provide an environment with seamless experimentation on a wide range of optimizations.


## Optimizations
Under the `src/baum_welch` directory, there are various different implementations of the Baum Welch algorithm. Each of those implementations contain a set of optimizations which aim to maximize the performance and minimize the running time of the algorithm. Below, we will analyze each of those implementations while also briefly specifying the exact optimizations that were employed:

### C++ implementations
* [`baum_baseline_impl.cc`](src/baum_welch/baum_baseline_impl.cc)
Naive and straightforward implementation using C++ STL, specifically `std::vector` was used extensively for all the data structures of the algorithm. This version is naive, though not deliberately underperforming. We eliminated obvious strided accesses and made sure to avoid obvious optimization blockers. Additionally, we avoided recursion through the dynamic programming paradigm.
* [`baum_baseline_opts.cc`](src/baum_welch/baum_baseline_opts.cc)
    * Numerical Optimizations - reduced the operations by performing mathematical analysis and finding out that we can eliminate the denominators in the expectation-maximization step of the algorithm. This can be trivially proven by induction.

### C-like implementations
* [`bw.cc`](src/baum_welch/baum_basic_opts.cc)
Naive and straightforward implementation using C-like pointers to represent all the arrays needed. This implementation was also implemented in a way to avoid obvious and major optimization blockers. This could be considered a pure C implementation.
* [`baum_basic_opts.cc`](src/baum_welch/baum_basic_opts.cc)
    * Reduced the operations by performing mathematical analysis and finding out that we can eliminate the denominators in the expectation-maximization step of the algorithm. This can be trivially proven by induction. (numerical optimization)
* [`bw_opts_v2.cc`](src/baum_welch/bw_opts_v2.cc)
    * Tweaked arrays' dimensions and access patterns to improve spatial and temporal locality (strided access elimination)
    * Eliminated redundant computations and used precomputation to avoid div operations inside loops
    * Eliminated auxiliary arrays and duplicated the computation due to heavy pollution of cache from 3D array
    * Minimized memory references in critical paths since we are memory bound (scalar replacement)
    * Merged loops and moved code based on execution ports availability (code motion)
    * Converted computations to as incremental and "weaker" as possible (strength reduction)
    * Eliminated `if` statement from inside the triple loop in expectation-maximization step by using an inverse map with `std::vector`. This optimization gave a huge bump to performance since we avoid a conditional execution which is input dependent and makes it impossible for branch predictor to achieve a decent hit rate.
* [`bw_loop_unroll.cc`](src/baum_welch/bw_loop_unroll.cc)
    * Optimally unrolled critical loops by using an autotuning infrastructure based on macros. Exploited ILP to the maximum.
* [`bw_opts_blocking.cc`](src/baum_welch/bw_opts_blocking.cc)
    * Added loop blocking techniques to improve the temporal and spatial locality of data accesses.
* [`bw_vectorized_1.cc`](src/baum_welch/bw_vectorized_1.cc)
    * Naive *vectorized* version by using Intel's vector intrinsics library and AVX2.
* [`bw_vectorized_opt.cc`](src/baum_welch/bw_vectorized_opt.cc)
    * Fully optimized *vectorized* version with heavy use of *FMAs*, unrolling and high execution port utilization.


## Building Procedure
### System Setup
In our project, we use CMake as our build configuration system. Particularly, the building of our
system requires:

* A C++17-enabled compiler. On Linux, gcc 8.3 should be sufficient.
* CMake 3.10 or higher
* On Linux and macOS, `make` utilities

On Ubuntu/Debian you can install the requirements with:

```
sudo apt-get install \
	  build-essential \
	  cmake
```

### Linux/Unix Build
Run from team32 directory after cloning the repo:
```
mkdir build && cd build
cmake ..
make
```
