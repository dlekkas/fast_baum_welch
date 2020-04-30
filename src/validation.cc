#include <iostream> //needed?
#include <assert.h>
#include <math.h>
#include "../include/validation.h"

//Error = Infinity Norm
double compute_error(double **transition_basic, double **emission_basic, double *init_dist_basic, double **transition_new,
  double **emission_new, double *init_dist_new, int M, int N) { //pointer parameters ok?

    int i, j;
    double diff;
    double max_diff= 0.0;

    for (i = 0; i < M; i++) {
        diff = std::abs(init_dist_new[i] - init_dist_basic[i]);
        if (diff > max_diff) {
          max_diff = diff;
        }
    }

    /*
    if (isnan(max_diff_pi)) {
      max_diff_pi = INFINITY;
    }*/

    for (i = 0; i < M; i++) {
        for (j = 0; j < M; j++) {
            diff = std::abs(transition_new[i][j] - transition_basic[i][j]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }
    }

    /*
    if (isnan(max_diff_transition)) {
      max_diff_transition = INFINITY;
    }*/

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            diff = std::abs(emission_new[i][j] - emission_basic[i][j]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }
    }

    /*
    if (isnan(max_diff_emission)) {
      max_diff_emission = INFINITY;
    }*/

    return max_diff;
  }
