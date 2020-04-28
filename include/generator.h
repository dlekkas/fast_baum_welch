/*
 * Hidden markov models are generative models, in which the joint distribution of observations
 * and hidden states, or equivalently both the prior distribution of hidden states (transition
 * probabilities) and conditional distribution of observations given states (emission probs),
 * is modeled. The algorithms operating on HMMs assume a prior uniform distribution over the
 * transition probabilities. However, it is also possible to create Hidden Markov Models with
 * other types of prior distributions. An obvious candidate, given the categorical distribution
 * of the transition probabilities, is the Dirichlet distribution, which is the conjugate prior
 * distribution of categorical distribution. Typically, a symmetric Dirichlet distribution is
 * chosen, reflecting ignorance about which states are inherently more likely than others.
 */

#ifndef GENERATOR_H
#define GENERATOR_H

#include <vector>

/*
 * Implementation based on Dirichlet Distribution from Wikipedia in the
 * subsection Gamma Distribution of section Random Number Generation:
 * https://en.wikipedia.org/wiki/Dirichlet_distribution#Random_number_generation
 */
std::vector<double> dirichlet_sample(int n_samples, const std::vector<double>& alpha, double beta = 1.0);

/* The normalization factor is expressed in terms of the multiplication
 * of gamma functions (see wiki link). For symmetric dirichlet distribution,
 * all elements making up the parameter vector alpha have the same value:
 * (https://en.wikipedia.org/wiki/Dirichlet_distribution). This single
 * alpha (i.e concentration parameter) controls the relative density or
 * sparseness of the resulting matrix
 */
std::vector<double> symmetric_dirichlet_sample(int n_samples, double concentrate = 1.0, double beta = 1.0);

std::vector<double> uniform_dist_sample(int n_samples);

std::vector<int> weighted_emission_sample(int n_samples, const std::vector<int>& weights);

std::vector<int> uniform_emission_sample(int n_samples, int M);


#endif
