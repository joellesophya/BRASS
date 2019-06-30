#ifndef COMPUTE_RESIDUALS_H
#define COMPUTE_RESIDUALS_H

#ifdef __cplusplus
extern "C" {
#endif

#define MISSVAL_FUNC_XI -999999

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "nrutil.h"
#include "brass.h"

typedef struct {
	int totpc;
	double **V;
	double *L;
	double *S;
	struct FAMILY *family;
	double *current_beta;
	int n_fam;
	int n_cov;
} param_func_xi;
typedef struct {
	param_func_xi parameters;
	double(*func_of_xi)(double, param_func_xi);
} function_xi;
double eqn3_lhs_minus_rhs(double xi, param_func_xi parameters);

double compute_gee_res(int method, struct FAMILY *family, struct DATA_STRUCT data_struct);

#ifdef __cplusplus
}
#endif

#endif 
