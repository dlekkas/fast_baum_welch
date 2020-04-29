#ifndef VALIDATION_H
#define VALIDATION_H

#define ERROR_BOUND 1e-2 // May Change it

double compute_error(double **transition_basic, double **emission_basic, double *init_dist_basic, double **transition_new,
  double **emission_new, double *init_dist_new, int M, int N);

//bool validation();

#endif
