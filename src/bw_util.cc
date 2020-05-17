#include <vector>
#include "../include/baum_welch.h"

using namespace std;



void copy_vec(const vector<vector<double>>& vec, double** arr) {
	for (size_t i = 0; i < vec.size(); i++) {
		copy(vec[i].begin(), vec[i].end(), arr[i]);
	}
}

void copy_vec(const vector<vector<double>>& vec, double* arr, int X) {
	for (size_t i = 0; i < vec.size(); i++) {
		copy(vec[i].begin(), vec[i].end(), arr + i*X);
	}
}

vector<double> convert_arr_to_vec(double* arr, int N) {
	vector<double> res(N);
	for (int i = 0; i < N; i++) {
		res[i] = arr[i];
	}
	return res;
}

vector<vector<double>> convert_arr_to_vec(double** arr, int X, int Y) {
	vector<vector<double>> res(X, vector<double>(Y));
	for (int i = 0; i < X; i++) {
		res[i] = convert_arr_to_vec(arr[i], Y);
	}
	return res;
}

vector<vector<double>> convert_arr_to_vec(double* arr, int X, int Y) {
	vector<vector<double>> res(X, vector<double>(Y));
	for (int i = 0; i < X; i++) {
		res[i] = convert_arr_to_vec(arr + i*Y, Y);
	}
	return res;
}



double* allocate_1d(int X) {
	return new double[X];
}


double** allocate_2d(int X, int Y) {
	double** ptr = new double*[X];
	for (auto i = 0; i < X; i++) {
		ptr[i] = allocate_1d(Y);
	}
	return ptr;
}

double*** allocate_3d(int X, int Y, int Z) {
	double*** ptr = new double**[X];
	for (auto i = 0; i < X; i++) {
		ptr[i] = allocate_2d(Y, Z);
	}
	return ptr;
}


double* allocate_alligned(int X, int Y) {
	double* ptr = (double *)aligned_alloc(32, X * Y * sizeof(double));
	return ptr;
}


void free_1d(double* arr) {
	if (arr != nullptr) {
		delete[] arr;
	}
}


void free_2d(double** arr, int X) {
	if (arr != nullptr) {
		for (auto i = 0; i < X; i++) {
			free_1d(arr[i]);
		}
	}
}


void free_3d(double*** arr, int X, int Y) {
	if (arr != nullptr) {
		for (auto i = 0; i < X; i++) {
			free_2d(arr[i], Y);
		}
	}
}



HMM BaumWelchC::GetHMM() {
	vector<vector<double>> tr = convert_arr_to_vec(A, M, M);
	vector<vector<double>> em = convert_arr_to_vec(B, N, M);
	vector<double> init_prob = convert_arr_to_vec(pi, M);
	HMM res(tr, em, init_prob);
	return res;
}

HMM BaumWelchCVect::GetHMM() {
	vector<vector<double>> tr = convert_arr_to_vec(A, M, M);
	vector<vector<double>> em = convert_arr_to_vec(B, N, M);
	vector<double> init_prob = convert_arr_to_vec(pi, M);
	HMM res(tr, em, init_prob);
	return res;
}


void BaumWelchC::Load(HMM& hmm, vector<int>& obs_seq) {
	M = hmm.transition.size();
	N = hmm.emission.size();
	T = obs_seq.size();

	A = allocate_2d(M, M);
	copy_vec(hmm.transition, A);

	B = allocate_2d(N, M);
	copy_vec(hmm.emission, B);

	pi = allocate_1d(M);
	copy(hmm.pi.begin(), hmm.pi.end(), pi);

	obs = new int[T];
	copy(obs_seq.begin(), obs_seq.end(), obs);

	fwd = allocate_2d(T, M);
	bwd = allocate_2d(T, M);
}

BaumWelchC::~BaumWelchC() {
	free_2d(A, M);
	free_2d(B, N);
	free_1d(pi);
	free_2d(fwd, T);
	free_2d(bwd, T);
	delete[] obs;
}


void BaumWelchCExtended::Load(HMM& hmm, vector<int>& obs_seq) {
	BaumWelchC::Load(hmm, obs_seq);
	gamma = allocate_2d(M, T);
	chsi = allocate_3d(M, M, T);
}

BaumWelchCExtended::~BaumWelchCExtended() {
	free_2d(gamma, M);
	free_3d(chsi, M, M);
}


void BaumWelchCVect::Load(HMM& hmm, vector<int>& obs_seq) {
	M = hmm.transition.size();
	N = hmm.emission.size();
	T = obs_seq.size();

	A = allocate_alligned(M, M);
	copy_vec(hmm.transition, A, M);

	B = allocate_alligned(N, M);
	copy_vec(hmm.emission, B, M);

	pi = allocate_alligned(M, 1);
	copy(hmm.pi.begin(), hmm.pi.end(), pi);

	obs = new int[T];
	copy(obs_seq.begin(), obs_seq.end(), obs);

	fwd = allocate_alligned(T, M);
	bwd = allocate_alligned(T, M);
}

BaumWelchCVect::~BaumWelchCVect() {

	free_1d(A);
	free_1d(B);
	free_1d(pi);
	free_1d(fwd);
	free_1d(bwd);
	delete[] obs;
}