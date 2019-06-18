#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

extern void getSeries(int dataPoints, double dt, int nComps, int nTypes, int* epsilon, int capacity, int* curState, int rng_seed, int* output);
