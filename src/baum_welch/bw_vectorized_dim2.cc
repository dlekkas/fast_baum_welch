/*
* Straightforward implementation of the Baum-Welch
* https://en.wikipedia.org/wiki/Baum–Welch_algorithm
*/
#include <iostream>
#include <assert.h>
#include <immintrin.h>

#include "../include/baum_welch.h"

using namespace std;


inline void transpose_square_matrix(double* matrix, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            swap(matrix[i*N + j], matrix[j*N + i]);
		}
	}
}

inline double hsum_double_avx(__m256d v) {
    __m128d vlow  = _mm256_castpd256_pd128(v);
    __m128d vhigh = _mm256_extractf128_pd(v, 1); // high 128
            vlow  = _mm_add_pd(vlow, vhigh);     // reduce down to 128
    __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
    return  _mm_cvtsd_f64(_mm_add_sd(vlow, high64));  // reduce to scalar
}

inline __m256d transpose_addition(__m256d a, __m256d b, __m256d c, __m256d d) {

    __m256d t1 = _mm256_unpacklo_pd(a, b);
    __m256d z1 = _mm256_unpackhi_pd(a, b);
    __m256d t2 = _mm256_unpacklo_pd(c, d);
    __m256d z2 = _mm256_unpackhi_pd(c, d);

    __m256d e1 = _mm256_add_pd(t1, z1);
    __m256d e2 = _mm256_add_pd(t2, z2);

    __m256d f1 = _mm256_permute4x64_pd(e1, 0b11101110);
    __m256d f2 = _mm256_permute4x64_pd(e2, 0b01000100);

    __m256d h1 = _mm256_add_pd(e1, f1);
    __m256d h2 = _mm256_add_pd(e2, f2);

    return _mm256_blend_pd(h1, h2, 0b1100);
}


inline __m256d transpose_addition_v2(__m256d a, __m256d b, __m256d c, __m256d d) {
	__m256d tmp0 = _mm256_hadd_pd(a, b);
    __m256d tmp1 = _mm256_hadd_pd(c, d);
    __m256d inter0 = _mm256_permute2f128_pd(tmp0, tmp1, 0x20);
    __m256d inter1 = _mm256_permute2f128_pd(tmp0, tmp1, 0x31);

	return _mm256_add_pd(inter0, inter1);
}





void forward_backward(double* forward, double* backward, int M, int N, int T,
		double* pi, double* A, double* B, int* observation_seq, double *sc_factors) {

	__m256d svec = _mm256_setzero_pd();
    for (int i = 0; i < M; i+=4) {
		__m256d b_vec = _mm256_load_pd(B + observation_seq[0]*N + i);
		__m256d pi_vec = _mm256_load_pd(pi + i);
		svec = _mm256_fmadd_pd(b_vec, pi_vec, svec);
		_mm256_store_pd(forward + i, _mm256_mul_pd(b_vec, pi_vec));
    }

    sc_factors[0] = 1.0 / hsum_double_avx(svec);


	// ops = M, mem = 2*M
	__m256d sc_vec = _mm256_set1_pd(hsum_double_avx(svec));
    for (int i = 0; i < M; i+=4) {
		__m256d fwd_vec = _mm256_load_pd(forward + i);
		_mm256_store_pd(forward + i, _mm256_mul_pd(fwd_vec, sc_vec));
    }


	transpose_square_matrix(A, M);

	// ops = 2*T*M^2 , mem = 2*T*M^2
    for (int t = 0; t < T-1; t++) {
		__m256d sum_vec = _mm256_setzero_pd();
	    for (int i = 0; i < M-3; i+=4) {
			__m256d b_vec1 = _mm256_load_pd(B + observation_seq[t+1]*M + i);
			__m256d accs1 = _mm256_setzero_pd(), accs2 = _mm256_setzero_pd();
			__m256d accs3 = _mm256_setzero_pd(), accs4 = _mm256_setzero_pd();

            for (int j = 0; j < M-3; j+=4) {
				__m256d fwd_vec1 = _mm256_load_pd(forward + t*M + j);

				__m256d A_vec1 = _mm256_load_pd(A + i*M + j);
				accs1 = _mm256_fmadd_pd(fwd_vec1, A_vec1, accs1);
				__m256d A_vec2 = _mm256_load_pd(A + (i+1)*M + j);
				accs2 = _mm256_fmadd_pd(fwd_vec1, A_vec2, accs2);
				__m256d A_vec3 = _mm256_load_pd(A + (i+2)*M + j);
				accs3 = _mm256_fmadd_pd(fwd_vec1, A_vec3, accs3);
				__m256d A_vec4 = _mm256_load_pd(A + (i+3)*M + j);
				accs4 = _mm256_fmadd_pd(fwd_vec1, A_vec4, accs4);
			}

			__m256d inter = transpose_addition(accs1, accs2, accs3, accs4);
			__m256d res1 = _mm256_mul_pd(b_vec1, inter);
			sum_vec = _mm256_add_pd(sum_vec, res1);
			_mm256_store_pd(forward + (t+1)*M + i, res1);

        }
        sc_factors[t+1] = 1.0 / hsum_double_avx(sum_vec);

		__m256d sc_vec = _mm256_set1_pd(sc_factors[t+1]);
        for (int i = 0; i < M; i+=8) {
			__m256d fwd_vec1 = _mm256_load_pd(forward + (t+1)*M + i);
			__m256d fwd_vec2 = _mm256_load_pd(forward + (t+1)*M + i+4);
			_mm256_store_pd(forward + (t+1)*M + i, _mm256_mul_pd(sc_vec, fwd_vec1));
			_mm256_store_pd(forward + (t+1)*M + i+4, _mm256_mul_pd(sc_vec, fwd_vec2));
		}
    }

	transpose_square_matrix(A, M);

	__m256d sum_vec = _mm256_set1_pd(sc_factors[T-1]);
    for (int i = 0; i < M; i+=4) {
		_mm256_store_pd(backward + (T-1)*M + i, sum_vec);
	}


    for (int t = T-2; t >= 0; t--) {
		__m256d sc_vec = _mm256_set1_pd(sc_factors[t]);
        for (int i = 0; i < M-3; i+=4) {
			__m256d acc1 = _mm256_setzero_pd(), acc2 = _mm256_setzero_pd();
			__m256d acc3 = _mm256_setzero_pd(), acc4 = _mm256_setzero_pd();
            for (int j = 0; j < M; j+=4) {
				__m256d bwd_vec = _mm256_load_pd(backward + (t+1)*M + j);
				__m256d b_vec= _mm256_load_pd(B + observation_seq[t+1]*M + j);
				__m256d inter = _mm256_mul_pd(bwd_vec, b_vec);

				__m256d A_vec1 = _mm256_load_pd(A + i*M + j);
				acc1 = _mm256_fmadd_pd(A_vec1, inter, acc1);
				__m256d A_vec2 = _mm256_load_pd(A + (i+1)*M + j);
				acc2 = _mm256_fmadd_pd(A_vec2, inter, acc2);
				__m256d A_vec3 = _mm256_load_pd(A + (i+2)*M + j);
				acc3 = _mm256_fmadd_pd(A_vec3, inter, acc3);
				__m256d A_vec4 = _mm256_load_pd(A + (i+3)*M + j);
				acc4 = _mm256_fmadd_pd(A_vec4, inter, acc4);
			}

			__m256d inter = transpose_addition(acc1, acc2, acc3, acc4);
			_mm256_store_pd(backward + t*M + i, _mm256_mul_pd(inter, sc_vec));
        }
    }

}


static void update_and_check(double* forward, double* backward, int M, int N, int T,
		double* pi, double* A, double* B, int* observation_seq, double *sc_factors,
		const vector<vector<int>>& obs_dict) {

	for (int t = 0; t < T; t+=4) {
		__m256d sc_vec = _mm256_load_pd(sc_factors+t);
		_mm256_store_pd(sc_factors + t, _mm256_div_pd(_mm256_set1_pd(1.0), sc_vec));
	}

	for (int z = 0; z < M; z+=4) {
		__m256d fwd_vec = _mm256_load_pd(forward + z);
		__m256d bwd_vec = _mm256_load_pd(backward + z);
		__m256d res = _mm256_mul_pd(fwd_vec, bwd_vec);
		_mm256_store_pd(pi + z, _mm256_mul_pd(res, _mm256_set1_pd(sc_factors[0])));
	}


	double temp_array[4];

	for (int i = 0; i < M; i+=4) {
		__m256d sum_vec = _mm256_setzero_pd();
		for (int t = 0; t < T-1; t++) {
            __m256d sc_vector = _mm256_set1_pd(sc_factors[t]);
            __m256d f_vec = _mm256_load_pd(forward + t*M + i);
            __m256d b_vec = _mm256_load_pd(backward + t*M + i);
			__m256d inter = _mm256_mul_pd(f_vec, b_vec);
			sum_vec = _mm256_fmadd_pd(inter, sc_vector, sum_vec);
		}
		sum_vec = _mm256_div_pd(_mm256_set1_pd(1.0) , sum_vec);
		_mm256_store_pd(temp_array, sum_vec);
		__m256d f1 = _mm256_set1_pd(temp_array[0]);
		__m256d f2 = _mm256_set1_pd(temp_array[1]);
		__m256d f3 = _mm256_set1_pd(temp_array[2]);
		__m256d f4 = _mm256_set1_pd(temp_array[3]);


		for (int j = 0; j < M-3; j+=4) {
			__m256d acc1 = _mm256_setzero_pd();
			__m256d acc2 = _mm256_setzero_pd();
			__m256d acc3 = _mm256_setzero_pd();
			__m256d acc4 = _mm256_setzero_pd();
			for (int t = 0; t < T-1; t++) {
				__m256d fwd_vec1 = _mm256_set1_pd(forward[t*M + i]);
				__m256d fwd_vec2 = _mm256_set1_pd(forward[t*M + i+1]);
				__m256d fwd_vec3 = _mm256_set1_pd(forward[t*M + i+2]);
				__m256d fwd_vec4 = _mm256_set1_pd(forward[t*M + i+3]);

				__m256d bwd_vec1 = _mm256_load_pd(backward + (t+1)*M + j);
				__m256d b_vec1 = _mm256_load_pd(B + observation_seq[t+1]*M + j);
				__m256d inter1 = _mm256_mul_pd(bwd_vec1, b_vec1);

				acc1 = _mm256_fmadd_pd(fwd_vec1, inter1, acc1);
				acc2 = _mm256_fmadd_pd(fwd_vec2, inter1, acc2);
				acc3 = _mm256_fmadd_pd(fwd_vec3, inter1, acc3);
				acc4 = _mm256_fmadd_pd(fwd_vec4, inter1, acc4);
			}
			acc1 = _mm256_mul_pd(acc1, f1);
			acc2 = _mm256_mul_pd(acc2, f2);
			acc3 = _mm256_mul_pd(acc3, f3);
			acc4 = _mm256_mul_pd(acc4, f4);
			_mm256_store_pd(A + i*M + j, _mm256_mul_pd(acc1, _mm256_load_pd(A + i*M + j)));
			_mm256_store_pd(A + (i+1)*M + j, _mm256_mul_pd(acc2, _mm256_load_pd(A + (i+1)*M + j)));
			_mm256_store_pd(A + (i+2)*M + j, _mm256_mul_pd(acc3, _mm256_load_pd(A + (i+2)*M + j)));
			_mm256_store_pd(A + (i+3)*M + j, _mm256_mul_pd(acc4, _mm256_load_pd(A + (i+3)*M + j)));
		}
	}



    for (int k = 0; k < M-3; k+=4) {
		__m256d sum1 = _mm256_setzero_pd();
        for (int t = 0; t < T; t++) {
			__m256d bwd_vec1 = _mm256_load_pd(backward + t*M + k);
			__m256d fwd_vec1 = _mm256_load_pd(forward + t*M + k);
			__m256d sc_vec1 = _mm256_set1_pd(sc_factors[t]);
			sum1 = _mm256_fmadd_pd(sc_vec1, _mm256_mul_pd(bwd_vec1, fwd_vec1), sum1);
		}

		for (int j = 0; j < N; j++) {
			__m256d acc1 = _mm256_setzero_pd();
			for (const auto& t: obs_dict[j]) {
				__m256d bwd_vec1 = _mm256_load_pd(backward + t*M + k);
				__m256d fwd_vec1 = _mm256_load_pd(forward + t*M + k);
				__m256d inter1 = _mm256_mul_pd(bwd_vec1, fwd_vec1);
				__m256d sc_vec1 = _mm256_set1_pd(sc_factors[t]);
				acc1 = _mm256_fmadd_pd(sc_vec1, inter1, acc1);
			}
			_mm256_store_pd(B + j*M + k, _mm256_div_pd(acc1, sum1));
		}
    }


}

void BaumWelchCVectDim2::operator()() {

    double *sc_factors = (double *)aligned_alloc(32, T * sizeof(double));

	vector<vector<int>> obs_dict(N);
	for (int t = 0; t < T; t++) {
		obs_dict[obs[t]].push_back(t);
	}

	for (int i = 0; i < MAX_ITERATIONS; i++) {
        forward_backward(fwd, bwd, M, N, T, pi, A, B, obs, sc_factors);
        update_and_check(fwd, bwd, M, N, T, pi, A, B, obs, sc_factors, obs_dict);
    }
}
