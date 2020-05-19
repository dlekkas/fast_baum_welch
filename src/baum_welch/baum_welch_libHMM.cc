#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../include/baum_welch.h"

// Code from: https://www.math.univ-toulouse.fr/~agarivie/Telecom/code/index.php

/*
sample from discrete distribution
in:   p = vector of probabilities,assumed to sum to 1
*/
int randm(double* p, int N){
	int res=0;
	double q=p[0];
	double u=(rand()+0.0)/RAND_MAX;
	while(u>q) q+=p[++res];
	return(res);
}

/*
sample a trajectory from a hidden markov chain
in:   nu = initial distribution as vector of size k
      Q = transition matrix of size k
      n = positive integer
out:  (x,y) = sample trajectory of size n of a HMM defined by (nu, Q, g):
      x = sample trajectory of size n of a Markov Chain with initial distribution nu and transition matrix Q
      y = observations such that the conditionnal distribution of y[k]
      given x[k] is g(x[k], :)
*/
void HMMsample(double* nu, double** Q, double** g, int* x, int* y, int n, int N, int M){
	int k;
	x[0] = randm(nu, N);
	y[0] = randm(g[x[0]], N);
	for(k=1; k<n; k++){
		x[k] = randm(Q[x[k-1]], N);
		y[k] = randm(g[x[k]] , M);
	}
}

/*
HMM filtering of an observation sequence, given hmm parameters
in:   y = vector of observations, assumed to be in range(g.shape[1])
      nu = initial distribution as vector of size k
      Q = transition matrix of size k x k
      g = emission matrix with k rows
out:  phi = filter: P(x[t]=x | y[0:t]=y[0:t]) for 0<=x<k and 0<=t<n
      c(t) = conditional likelihood: P(Y[t] = y[t]| Y[0:t-1]=y[0:t-1])
*/
void HMMfilter(int* y, double* nu, double** Q, double** g, double** phi, double* c, int n, int N, int M){
	int i,j,t;
	double z[N];
	c[0]=0;
	for(j=0; j<N; j++){
		z[j] = nu[j]*g[j][y[0]];
		c[0] += z[j];
	}
	for(j=0; j<N; j++) phi[0][j] = z[j]/c[0];
	for(t=1; t<n; t++){
		c[t]=0;
		for(j=0; j<N; j++){
			z[j]=0;
			for(i=0; i<N; i++) z[j]+=phi[t-1][i]*Q[i][j]*g[j][y[t]];
			c[t] += z[j];
		}
		for(j=0; j<N; j++) phi[t][j] = z[j]/c[t];
	}
}

/*
HMM filtering of an observation sequence, given hmm parameters
in:   y = vector of observations, assumed to be in range(Q.shape[0])
      Q = transition matrix of size k x k
      g = emission matrix with k rows
      c = conditional likelihoods, computed by HMMfilter
out:  beta = smoothing factors: P(y[t+1:n]=y[t+1:n] | X[t]=x) / P(Y[t+1:n]=y[t+1:n] | Y[0:t]=y[1:t]) for 0<=x<k and 1<=t<n
      permits to compute the posterior distribution of the hidden states
      P(X[t]=x | Y[0:n]=y[0:n])  as post = phi .* beta
*/
void HMMsmoother(int* y, double** Q, double** g, double* c, double** beta, int n, int N, int M){
	int i,j,t;
	for(j=0; j<N; j++) beta[n-1][j] = 1;
	for(t=n-2; t>=0; t--){
		for(i=0; i<N; i++){
			double z=0;
			for(j=0; j<N; j++)
				z += Q[i][j] * g[j][y[t+1]] * beta[t+1][j];
			beta[t][i] = z / c[t+1];
		}
	}
}

/*
utility functions: sample random transition and emission kernels
*/
void randomTransitionKernel(double** K, int N){
	int i,j;
	double s;
	for(i=0; i<N; i++){
		s=0;
		for(j=0; j<N; j++)
			s+=(K[i][j] = (rand()+0.0)/RAND_MAX);
		for(j=0; j<N; j++)
			K[i][j] /= s;
	}
}

void randomEmissionKernel(double** K, int N, int M){
	int i,j;
	double s;
	for(i=0; i<N; i++){
		s=0;
		for(j=0; j<M; j++)
			s+=(K[i][j] = (rand()+0.0)/RAND_MAX);
		for(j=0; j<M; j++)
			K[i][j] /= s;
	}
}

/*
compute maximum likehood estimate using Expectation-Maximization
iterations
in:   y = vector of observations
      nu = initial distribution of the hidden chain
      tol = tolerance for the stopping criterion
      maxIt = maximal number of iterations
out:  Q = estimate of the transition matrix of the hidden markov process
      g = estimated probabilities of transition: g(x,y) = estimate of P(Y=y | X=x) for 0<=x<k
      l = log-likelihood of y for parameters Q and g
*/
void BaumWelchLibHMM::operator()(){
	// returns best log-likelihood found
	// static const int maxIt = 10;
	// static const double tol = 1e-4;
	int i, j, it, t;
	double z, l=0;
	randomTransitionKernel(Q, N);
	randomEmissionKernel(g, N, M);
	double c[n];
	double A[N][M];
	double s[N];
	double B[N][N];
	double s2[N];
	for(it=0; (it<MAX_ITERATIONS); it++){
		//printf("Iteration: %d\n", it);
		HMMfilter(y, nu, Q, g, phi, c, n, N, M);
		HMMsmoother(y, Q, g, c, beta, n, N, M);
		for(i=0; i<N; i++){
			s[i]=0; s2[i]=0;
			for(j=0; j<M; j++) A[i][j] = 0;
			for(j=0; j<N; j++) B[i][j] = 0;
		}
		for(t=0; t<n; t++)
			for(i=0; i<N; i++){
				z=phi[t][i]*beta[t][i];
				A[i][y[t]] += z;
				s[i] += z;
			}
		for(i=0; i<N; i++)
			for(j=0; j<M; j++){
				g[i][j] = A[i][j] / s[i];
			}
		for(t=1; t<n; t++){
			for(i=0; i<N; i++)
				for(j=0; j<N; j++){
					z=phi[t-1][i]*Q[i][j]*g[j][y[t]]*beta[t][j]/c[t];
					B[i][j] += z;
					s2[i] += z;
				}
		}
		for(i=0; i<N; i++)
			for(j=0; j<N; j++){
				Q[i][j] = B[i][j] / s2[i];
			}
	}
}
