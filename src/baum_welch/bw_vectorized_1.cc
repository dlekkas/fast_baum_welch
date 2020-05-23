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

inline void forward_backward(double* forward, double* backward, int M, int N, int T,
		double* pi, double* A, double* B, int* observation_seq, double *sc_factors) {

    double sum = 0.0, acc = 0.0;

    __m256d B_vect, p, sum_vector, temp, f, b, A_vect, sc_vector;

	// ops = 2*M , mem = 4*M
    for (int i = 0; i < M; i+=4) {

        double* B0 = B + observation_seq[0]*N;
        B_vect = _mm256_load_pd(B0 + i);
        p = _mm256_load_pd(pi + i);
        temp = _mm256_mul_pd(p, B_vect);
        _mm256_store_pd(forward + i, temp);

        for (int j=0; j<4; j++)
            sum += forward[i+j];
    }
    sc_factors[0] = 1.0 / sum;

	// ops = M, mem = 2*M
    sum_vector = _mm256_set1_pd(sum);
    for (int i = 0; i < M; i+=4) {
        f = _mm256_load_pd(forward + i);
        temp = _mm256_mul_pd(f, sum_vector);
        _mm256_store_pd(forward + i, temp);
    }

    double temp_array[4];
	// ops = 2*T*M^2 , mem = 2*T*M^2
    for (int t = 1; t < T; t++) {
        sum = 0.0;
	    for (int i = 0; i < M; i++) {
            acc = 0.0;
            for (int j = 0; j < M; j+=4) {
                f = _mm256_load_pd(forward + (t-1)*M + j);
                A_vect = _mm256_set_pd(A[(j+3)*M + i], A[(j+2)*M + i], A[(j+1)*M + i], A[j*M + i]);
                temp = _mm256_mul_pd(f, A_vect);
                _mm256_store_pd(temp_array, temp);
                for (int k=0; k<4; k++)
                    acc += temp_array[k];
			}

            //acc = _mm256_cvtsd_f64(acc_v);
            forward[t*M + i] = B[observation_seq[t]*M + i] * acc;
            sum += forward[t*M + i];

        }

        sc_factors[t] = 1.0 / sum;
        sum_vector = _mm256_set1_pd(sc_factors[t]);
        for (int i = 0; i < M; i+=4) {
            f = _mm256_load_pd(forward + t*M + i);
            temp = _mm256_mul_pd(f, sum_vector);
            _mm256_store_pd(forward + t*M + i, temp);
		}

    }

    sc_vector = _mm256_set1_pd(sc_factors[T-1]);
    for (int i = 0; i < M; i+=4) {
        _mm256_store_pd(backward + (T-1)*M + i, sc_vector);
	}

	// ops = 3*T*M^2, mem = 3*T*M^2
    for (int t = T-2; t >= 0; t--) {
        for (int i = 0; i < M; i++) {
            sum = 0.0;
            for (int j = 0; j < M; j+=4) {
                b = _mm256_load_pd(backward + (t+1)*M + j);
                A_vect = _mm256_load_pd(A + i*M + j);
                B_vect = _mm256_load_pd(B + observation_seq[t+1]*M + j);

                temp = _mm256_mul_pd(b, A_vect);
                temp = _mm256_mul_pd(temp, B_vect);
                _mm256_store_pd(temp_array, temp);
                for (int k=0; k<4; k++)
                    sum += temp_array[k];
			}

            backward[t*M + i] = sum*sc_factors[t];
        }
    }

}


inline void update_and_check(double* forward, double* backward, int M, int N, int T,
		double* pi, double* A, double* B, int* observation_seq, double *sc_factors, vector<vector<int>>& obs_dict) {

    __m256d B_vect, sum_vector, temp, f, b, acc_v, sc_vector, occurences, acc1;

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
            temp = _mm256_mul_pd(temp, sc_vector);
            acc_v = _mm256_add_pd(acc_v, temp);
		}

		for (int j = 0; j < M; j++) {
            sum_vector = _mm256_set1_pd(0.0);
			for (int t = 0; t < T-1; t++) {
                B_vect = _mm256_set1_pd(backward[(t+1)*M + j] * B[(observation_seq[t+1])*M + j]);
                f = _mm256_load_pd(forward + t*M + i);
                temp = _mm256_mul_pd(B_vect, f);
                sum_vector = _mm256_add_pd(sum_vector, temp);
            }
            temp =  _mm256_div_pd(sum_vector, acc_v);
            _mm256_store_pd(temp_array, temp);

            A[i*M + j] *= temp[0];
            A[(i+1)*M + j] *= temp[1];
            A[(i+2)*M + j] *= temp[2];
            A[(i+3)*M + j] *= temp[3];

		}

    }
    for (int i = 0; i < M; i+=4) {

        sum_vector = _mm256_set1_pd(0.0);

        for (int t = 0; t < T; t++) {
            sc_vector = _mm256_set1_pd(1.0/sc_factors[t]);
            f = _mm256_load_pd(forward + t*M + i);
            b = _mm256_load_pd(backward + t*M + i);
            temp = _mm256_mul_pd(f, b);
            temp = _mm256_mul_pd(temp, sc_vector);
            sum_vector = _mm256_add_pd(temp, sum_vector);
		}


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

void BaumWelchCVectBasic::operator()() {

    vector<vector<int>> obs_dict(N);
	for (int t = 0; t < T; t++) {
		obs_dict[obs[t]].push_back(t);
	}

    double *sc_factors = (double *)aligned_alloc(32, T * sizeof(double));

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