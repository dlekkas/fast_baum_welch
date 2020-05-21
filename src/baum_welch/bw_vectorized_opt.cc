/*
* Straightforward implementation of the Baum-Welch
* https://en.wikipedia.org/wiki/Baum–Welch_algorithm
*/
#include <iostream>
#include <assert.h>
#include <immintrin.h>

#include "../include/baum_welch.h"

// forward, backward: T*M
// A: M*M
// B: N*M
// pi : M

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

inline __m256d add_4_arrays(__m256d a, __m256d b, __m256d c, __m256d d) {
      __m256d t1, t2, z1, z2, e1, e2, f1, f2, h1, h2;


    // first we create arrays: [A[1], B[1], A[3], B[3]], [A[2], B[2], A[4], B[4]], [C[1], D[1], C[3], D[3]], [C[2], D[2], C[4], D[4]]
    t1 = _mm256_unpacklo_pd(a, b);
    z1 = _mm256_unpackhi_pd(a, b);
    t2 = _mm256_unpacklo_pd(c, d);
    z2 = _mm256_unpackhi_pd(c, d);

    // now we create: [A[1]+A[2], B[1]+B[2], A[3]+A[4], B[3]+B[4]], [C[1]+C[2], D[1]+D[2], C[3]+C[4], D[3]+D[4]]
    e1 = _mm256_add_pd(t1, z1);
    e2 = _mm256_add_pd(t2, z2);

    // now we create: [A[3]+A[4], B[3]+B[4], A[3]+A[4], B[3]+B[4]], [C[3]+C[4], D[3]+D[4], C[3]+C[4], D[3]+D[4]]
    f1 = _mm256_permute4x64_pd(e1, 0b11101110);
    f2 = _mm256_permute4x64_pd(e2, 0b01000100);

    // now we create [A[1]+A[2]+A[3]+A[4], B[1]+B[2]+B[3]+B[4], .... ]
    h1 = _mm256_add_pd(e1, f1);

    // now we create [....., C[1]+C[2]+C[3]+C[4], D[1]+D[2]+D[3]+D[4]]
    h2 = _mm256_add_pd(e2, f2);

    // [A[1]+A[2]+A[3]+A[4], B[1]+B[2]+B[3]+B[4], C[1]+C[2]+C[3]+C[4], D[1]+D[2]+D[3]+D[4]]
    return _mm256_blend_pd(h1, h2, 0b1100);

}


inline __m256d transpose_addition(__m256d a, __m256d b, __m256d c, __m256d d) {
	__m256d tmp0 = _mm256_hadd_pd(a, b);
    __m256d tmp1 = _mm256_hadd_pd(c, d);
    __m256d inter0 = _mm256_permute2f128_pd(tmp0, tmp1, 0x20);
    __m256d inter1 = _mm256_permute2f128_pd(tmp0, tmp1, 0x31);

	return _mm256_add_pd(inter0, inter1);
}




inline void forward_backward(double* forward, double* backward, int M, int N, int T,
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



inline void update_and_check(double* forward, double* backward, int M, int N, int T,
		double* pi, double* A, double* B, int* observation_seq, double *sc_factors, vector<vector<int>>& obs_dict, double* accs) {

    __m256d B_vect, sum_vector, temp, f, b, acc_v, sc_vector;
    __m256d acc1;
    __m256d sum_vector1, sum_vector2, sum_vector3, sum_vector4;
    __m256d f1, f2, f3, f4;
    __m256d A1, A2, A3, A4;
    __m256d temp1, temp2, temp3, temp4;


    double temp_array[4];

	// ops = 3*T*M^2 , mem = 4*T*M^2
    sc_vector = _mm256_set1_pd(1.0/sc_factors[0]);
	for (int i = 0; i < M; i+=4) {

        f = _mm256_load_pd(forward + i);
        b = _mm256_load_pd(backward + i);

        temp = _mm256_mul_pd(f, b);
        temp = _mm256_mul_pd(temp, sc_vector);
        _mm256_store_pd(pi+i, temp);
    }

    for (int i=0; i<M; i+=4) {

        acc_v = _mm256_load_pd(pi + i);

        for (int t = 1; t < T-1; t++) {
            sc_vector = _mm256_set1_pd(1.0/sc_factors[t]);
            f = _mm256_load_pd(forward + t*M + i);
            b = _mm256_load_pd(backward + t*M + i);

            temp = _mm256_mul_pd(f, b);
            acc_v = _mm256_fmadd_pd(temp, sc_vector, acc_v);
             _mm256_store_pd(accs+i, acc_v);

		}

        for (int j = 0; j < M; j+=4) {
            sum_vector1 = _mm256_set1_pd(0.0);
            sum_vector2 = _mm256_set1_pd(0.0);
            sum_vector3 = _mm256_set1_pd(0.0);
            sum_vector4 = _mm256_set1_pd(0.0);

			for (int t = 0; t < T-1; t++) {
                b = _mm256_load_pd(backward + (t+1)*M + j);
                B_vect = _mm256_load_pd(B + observation_seq[t+1]*M + j);
                temp = _mm256_mul_pd(b, B_vect);

                f1 = _mm256_set1_pd(forward[t*M + i]);
                f2 = _mm256_set1_pd(forward[t*M + i+1]);
                f3 = _mm256_set1_pd(forward[t*M + i+2]);
                f4 = _mm256_set1_pd(forward[t*M + i+3]);

                sum_vector1 = _mm256_fmadd_pd(temp, f1, sum_vector1);
                sum_vector2 = _mm256_fmadd_pd(temp, f2, sum_vector2);
                sum_vector3 = _mm256_fmadd_pd(temp, f3, sum_vector3);
                sum_vector4 = _mm256_fmadd_pd(temp, f4, sum_vector4);
            }

            _mm256_store_pd(temp_array, acc_v);
            f1 = _mm256_set1_pd(1.0/temp_array[0]);
            f2 = _mm256_set1_pd(1.0/temp_array[1]);
            f3 = _mm256_set1_pd(1.0/temp_array[2]);
            f4 = _mm256_set1_pd(1.0/temp_array[3]);

            sum_vector1 = _mm256_mul_pd(sum_vector1, f1);
            sum_vector2 = _mm256_mul_pd(sum_vector2, f2);
            sum_vector3 = _mm256_mul_pd(sum_vector3, f3);
            sum_vector4 = _mm256_mul_pd(sum_vector4, f4);

            A1 = _mm256_load_pd(A + i*M + j);
            A2 = _mm256_load_pd(A + (i+1)*M + j);
            A3 = _mm256_load_pd(A + (i+2)*M + j);
            A4 = _mm256_load_pd(A + (i+3)*M + j);

            temp1 = _mm256_mul_pd(A1, sum_vector1);
            temp2 = _mm256_mul_pd(A2, sum_vector2);
            temp3 = _mm256_mul_pd(A3, sum_vector3);
            temp4 = _mm256_mul_pd(A4, sum_vector4);

            _mm256_store_pd(A + i*M + j, temp1);
            _mm256_store_pd(A + (i+1)*M + j, temp2);
            _mm256_store_pd(A + (i+2)*M + j, temp3);
            _mm256_store_pd(A + (i+3)*M + j, temp4);


        }
    }

	// ops = 3*T*M*N, mem = 3*T*M^2
    for (int i = 0; i < M; i+=4) {

        acc_v = _mm256_load_pd(accs+i);
        sc_vector = _mm256_set1_pd(1.0/sc_factors[T-1]);
        f = _mm256_load_pd(forward + (T-1)*M + i);
        b = _mm256_load_pd(backward + (T-1)*M + i);

        temp = _mm256_mul_pd(f, b);
        sum_vector = _mm256_fmadd_pd(temp, sc_vector, acc_v);

        for (int j=0; j < N; j++) {
            acc1 = _mm256_set1_pd(0.0);
			for (const auto& t: obs_dict[j]) {
                sc_vector = _mm256_set1_pd(1.0/sc_factors[t]);
				f = _mm256_load_pd(forward + t*M + i);
                b = _mm256_load_pd(backward + t*M + i);
                temp = _mm256_mul_pd(f, b);
                acc1 = _mm256_fmadd_pd(temp, sc_vector, acc1);
			}

			temp = _mm256_div_pd(acc1, sum_vector);
            _mm256_store_pd(B + j*M + i, temp);

        }
    }



}

void BaumWelchCVectOpt::operator()() {

    /* Significant optimization to get rid of conditional addition with
	 * extremely bad branch prediction patterns (highly unlikely to
	 * predict its outcome */
	vector<vector<int>> obs_dict(N);
	for (int t = 0; t < T; t++) {
		obs_dict[obs[t]].push_back(t);
	}

    double *sc_factors = (double *)aligned_alloc(32, T * sizeof(double));
    double *accs = (double *)aligned_alloc(32, M * sizeof(double));


	for (int i = 0; i < MAX_ITERATIONS; i++) {
        forward_backward(fwd, bwd, M, N, T, pi, A, B, obs, sc_factors);
        update_and_check(fwd, bwd, M, N, T, pi, A, B, obs, sc_factors, obs_dict, accs);
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
