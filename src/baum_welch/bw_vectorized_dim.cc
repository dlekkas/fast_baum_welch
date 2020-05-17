/*
* Straightforward implementation of the Baum-Welch
* https://en.wikipedia.org/wiki/Baum–Welch_algorithm
*/
#include <iostream>
#include <assert.h>
#include <immintrin.h>
//#include <x86intrin.h>

#include "../include/baum_welch.h"

// forward, backward: T*M
// A: M*M
// B: N*M
// pi : M

using namespace std;


inline double hsum_double_avx(__m256d v) {
    __m128d vlow  = _mm256_castpd256_pd128(v);
    __m128d vhigh = _mm256_extractf128_pd(v, 1); // high 128
            vlow  = _mm_add_pd(vlow, vhigh);     // reduce down to 128
    __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
    return  _mm_cvtsd_f64(_mm_add_sd(vlow, high64));  // reduce to scalar
}


static void forward_backward(double* forward, double* backward, int M, int N, int T,
		double* pi, double* A, double* B, int* observation_seq, double *sc_factors) {

	// ops = 2*M , mem = 4*M
	double sum = 0.0;
    for (int i = 0; i < M; i++) {
		forward[i] = pi[i] * B[observation_seq[0]*N + i];
		sum += forward[i];
    }
    sc_factors[0] = 1.0 / sum;

	// ops = M, mem = 2*M
    for (int i = 0; i < M; i++) {
		forward[i] *= sum;
    }


	// ops = 2*T*M^2 , mem = 2*T*M^2
    for (int t = 0; t < T-1; t++) {
        sum = 0.0;
	    for (int i = 0; i < M-7; i+=8) {
			__m256d accs1 = _mm256_setzero_pd(), accs2 = _mm256_setzero_pd();
			__m256d accs3 = _mm256_setzero_pd(), accs4 = _mm256_setzero_pd();
            for (int j = 0; j < M-1; j+=2) {
				__m256d A_vec1 = _mm256_load_pd(A + j*M + i);
				__m256d A_vec2 = _mm256_load_pd(A + j*M + i+4);

				__m256d fwd_vec1 = _mm256_set1_pd(forward[t*M + j]);
				accs1 = _mm256_fmadd_pd(fwd_vec1, A_vec1, accs1);
				__m256d fwd_vec2 = _mm256_set1_pd(forward[t*M + j+1]);
				accs2 = _mm256_fmadd_pd(fwd_vec1, A_vec2, accs2);
				__m256d A_vec3 = _mm256_load_pd(A + (j+1)*M + i);
				accs3 = _mm256_fmadd_pd(fwd_vec2, A_vec3, accs3);
				__m256d A_vec4 = _mm256_load_pd(A + (j+1)*M + i+4);
				accs4 = _mm256_fmadd_pd(fwd_vec2, A_vec4, accs4);
			}
			__m256d b_vec1 = _mm256_load_pd(B + observation_seq[t+1]*M + i);
			__m256d b_vec2 = _mm256_load_pd(B + observation_seq[t+1]*M + i+4);

			__m256d res1 = _mm256_mul_pd(b_vec1, _mm256_add_pd(accs1, accs3));
			_mm256_store_pd(forward + (t+1)*M + i, res1);
			__m256d res2 = _mm256_mul_pd(b_vec2, _mm256_add_pd(accs2, accs4));
			_mm256_store_pd(forward + (t+1)*M + i+4, res2);

			sum += hsum_double_avx(_mm256_add_pd(res1, res2));
        }
        sc_factors[t+1] = 1.0 / sum;
        for (int i = 0; i < M; i++) {
			forward[(t+1)*M + i] = forward[(t+1)*M + i] / sum;
		}
    }



    for (int i = 0; i < M; i++) {
		backward[(T-1)*M + i] = sc_factors[T-1];
	}

	// ops = 3*T*M^2, mem = 3*T*M^2
    for (int t = T-2; t >= 0; t--) {
        for (int i = 0; i < M; i++) {
            sum = 0.0;
            for (int j = 0; j < M; j++) {
				sum += A[i*M + j] * backward[(t+1)*M + j] * B[observation_seq[t+1]*M + j];
			}

            backward[t*M + i] = sum*sc_factors[t];
        }
    }

}


static void update_and_check(double* forward, double* backward, int M, int N, int T,
		double* pi, double* A, double* B, int* observation_seq, double *sc_factors,
		const vector<vector<int>>& obs_dict) {

	for (int t = 0; t < T; t++) {
		sc_factors[t] = 1.0 / sc_factors[t];
	}

	for (int z = 0; z < M; z++) {
		pi[z] = forward[z] * backward[z] * sc_factors[0];
	}



	// ops = 3*T*M^2 , mem = 4*T*M^2

	for (int i = 0; i < M; i++) {
		double acc = 0.0;
		for (int t = 0; t < T-1; t++) {
			acc += (forward[t*M + i] * backward[t*M + i]) * sc_factors[t];
		}
		for (int j = 0; j < M; j++) {
			double sum = 0.0;
			for (int t = 0; t < T-1; t++) {
				sum += backward[(t+1)*M + j] * B[observation_seq[t+1]*M + j] * forward[t*M + i];
			}
			A[i*M + j] *= sum / acc;
		}
	}



	// ops <= 4*T*M*N, mem <= 3*T*M^2
    for (int k = 0; k < M; k++) {
        double sum1 = 0.0;
        for (int t = 0; t < T; t++) {
            sum1 += forward[t*M + k] * backward[t*M + k] * sc_factors[t];
		}

		for (int j = 0; j < N; j++) {
			double acc = 0.0;
			for (const auto& t: obs_dict[j]) {
				acc += (forward[t*M + k] * backward[t*M + k]) * sc_factors[t];
			}
			B[j*M + k] = acc / sum1;
		}
    }


}

void BaumWelchCVectDim::operator()() {

    double *sc_factors = (double *)aligned_alloc(32, T * sizeof(double));

	vector<vector<int>> obs_dict(N);
	for (int t = 0; t < T; t++) {
		obs_dict[obs[t]].push_back(t);
	}


	for (int i = 0; i < MAX_ITERATIONS; i++) {
        forward_backward(fwd, bwd, M, N, T, pi, A, B, obs, sc_factors);
        update_and_check(fwd, bwd, M, N, T, pi, A, B, obs, sc_factors, obs_dict);
    }

	#ifdef DEBUG
		using namespace std;

		cout << endl << "--------------------------  DEBUG START -------------------------- " << endl;
		cout << "Matrix A:" << endl;
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				cout << A[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl << "Matrix B:" << endl;
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < N; j++) {
				cout << B[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl << "Pi vector:" << endl;
		for (int i = 0; i < M; i++) {
			cout << pi[i] << " ";
		}
		cout << endl << endl;
		cout << endl << "-------------------------- DEBUG END -------------------------- " << endl;
	#endif

}
