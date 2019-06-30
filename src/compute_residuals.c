#include "compute_residuals.h"

double compute_gee_res(int method, struct FAMILY *family, struct DATA_STRUCT data_struct){

    int i, j, fam = 1;
    int n_total = data_struct.n_total;
    int n_cov = data_struct.n_cov;
    int n_fam = data_struct.n_fam;
    double *yaffec, *Dinv, *optim_betahat, *lambda , **u, **Cov, optim_xi = -1.0;

    // Arrays to fit CARAT null model and sample trait replicates
    if(method == 0) {
        (family + fam)->trait = dvector(1, n_total);
        (family + fam)->beta_hat = dvector(1, n_cov);
        (family + fam)->Dinv = dvector(1, n_total);
        (family + fam)->mu_hat_QL = dvector(1, n_total);
        (family + fam)->pre_mult_mat = dmatrix(1, n_total, 1, n_total - n_cov);
        (family + fam)->delta = dvector(1, n_total - n_cov);
    }

    Cov = (family+fam)->cov;
    lambda = (family + fam)->eigvalue;
    u = (family + fam)->eigvector;
    optim_betahat = (family + fam)->beta_hat;
    Dinv = (family + fam)->Dinv;

    switch(method){
        case 0: {
            yaffec = (family + fam)->trait_obs;
            break;
        }
        case 1: case 2:{
            yaffec = (family + fam)->trait_rep;
            break;
        }
        default:
            printf("Wrong method (%d) indicated!\n", method);
            exit(1);
    }

    // Current response analyzed
    for (j = 1; j <= n_total; j++)
        (family+fam)->trait[j] = yaffec[j];

    // Start model fitting
    //clock_t time1 = clock();
    compute_logistic_coeff(yaffec, Cov, n_total, n_cov, optim_betahat);

    estimate_gee_parameters(method, n_total, u, lambda, Dinv, family, n_cov, n_fam, &optim_xi, optim_betahat);
    //clock_t time2 = clock();
    //printf("Time for fitting CARAT null model: %lf\n", (double)(time2 - time1)/ CLOCKS_PER_SEC);

    return optim_xi;
}

double eqn3_lhs_minus_rhs(double xi, param_func_xi parameters) {
    return target_func_xi(xi, parameters.totpc, parameters.V, parameters.L, parameters.S, parameters.family, parameters.current_beta, parameters.n_fam, parameters.n_cov);
}
