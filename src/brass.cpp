#include <Eigen>
#include <algorithm> // For shuffling
#include <string> 
#include <stdio.h>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "iostream"
#include "brass.h"

using Eigen::MatrixXd;
using namespace Eigen;
using namespace std;

typedef Map<Matrix<double, Dynamic, Dynamic, RowMajor> > MapRowXd;
typedef Map<Matrix<int, Dynamic, Dynamic, RowMajor> > MapRowXi;
typedef Map<Matrix<double, Dynamic, Dynamic, ColMajor> > MapColXd;
typedef Map<ArrayXXd> MapArXd;

#ifdef __cplusplus
extern "C"
{
#endif

    // Lapack eigen-decomposition
    void dsyevr_( char* jobz, char* range, char* uplo, int* n, double* a, int* lda, double* vl, double* vu, int* il, int* iu, double* abstol, int* M,
            double* w, double* z, int* ldz, int* isuppz, double* work, int* lwork, int* iwork, int* liwork, int* info );

#include "nrutil.h"
#include "nrutil.c"
#include "brent.h"
#include "hashes.h"
#include "datacheck.h"
#include "read.h"
#include "compute_residuals.h"
#include "permute_sim.h"

#ifdef __cplusplus
}
#endif

int main(int argc, char **argv){

    printf("\n"
            "---------------------------------------------------------\n"
            "|                         BRASS                         |\n"
            "|                                                       |\n"
            "|              Binary trait Resampling method           |\n"
            "|              Adjusting for Sample Structure           |\n"
            "|                                                       |\n"
            "|             Version 1.0 - December 07, 2018           |\n"
            "|                                                       |\n"
            "|                  Copyright(C) 2018-2020               |\n"
            "|     Joelle Mbatchou, Mark Abney and Mary Sara McPeek  |\n"
            "|                                                       |\n"
            "|                        Homepage:                      |\n"
            "|     http://galton.uchicago.edu/~mcpeek/software/BRASS |\n"
            "---------------------------------------------------------\n\n"
          );


    // Files in
    char pheno_filename[MAXFILELEN]= "pheno.txt", eig_filename[MAXFILELEN], kin_filename[MAXFILELEN], *end = NULL;
    // Files out
    char out_prefix[MAXFILELEN], param_filename[MAXFILELEN], yrep_filename[MAXFILELEN], outeig_filename[MAXFILELEN];
    // For options
    int pfile = 0, efile = 0, kfile = 0, ofile = 0;
    int lapa = 1, print_eig = 0, nperm_in = 0, n_perm, seed_in = 0;
    long nperm_long, seed_rand;

    for (int arg = 1; arg < argc; arg++) {
        if (strcmp(argv[arg], "-p") == 0 || strcmp(argv[arg], "--pheno") == 0){
            if (argv[arg + 1] == NULL || argv[arg + 1][0] == '-') {
                continue;
            }
            arg++;
            pheno_filename[MAXFILELEN-1]!='\0';
            strncpy(pheno_filename, argv[arg], MAXFILELEN);
            if(pheno_filename[MAXFILELEN-1]!='\0'){
                printf("ERROR: Input file names cannot exceed %d characters. Please double check the phenotype and covariate data file!\n", MAXFILELEN-1);
                exit(1);
            }
            pfile = 1;
        } else if (strcmp(argv[arg], "-e") == 0 || strcmp(argv[arg], "--eigen") == 0){
            if (argv[arg + 1] == NULL || argv[arg + 1][0] == '-') {
                continue;
            }
            arg++;
            eig_filename[MAXFILELEN-1]='\0';
            strncpy(eig_filename, argv[arg], MAXFILELEN);
            if(eig_filename[MAXFILELEN-1]!='\0'){
                printf("ERROR: Input file names cannot exceed %d characters. Please double check the eigen-decomposition file!\n", MAXFILELEN-1);
                exit(1);
            }
            efile = 1;
        } else if (strcmp(argv[arg], "-E") == 0){
            lapa = 0;
        } else if (strcmp(argv[arg], "-k") == 0 || strcmp(argv[arg], "--kin") == 0){
            if (argv[arg + 1] == NULL || argv[arg + 1][0] == '-') {
                continue;
            }
            arg++;
            kin_filename[MAXFILELEN-1]='\0';
            strncpy(kin_filename, argv[arg], MAXFILELEN);
            if(kin_filename[MAXFILELEN-1]!='\0'){
                printf("ERROR: Input file names cannot exceed %d characters. Please double check the GRM file!\n", MAXFILELEN-1);
                exit(1);
            }
            kfile = 1;
        } else if (strcmp(argv[arg], "-o") == 0 || strcmp(argv[arg], "--out") == 0){
            if (argv[arg + 1] == NULL || argv[arg + 1][0] == '-') {
                continue;
            }
            arg++;
            strncpy(out_prefix, argv[arg], MAXFILELEN);
            ofile = 1;
            if (strlen(out_prefix)>MAXFILELEN-10) {
                printf("ERROR: The prefix for the output files cannot be longer than %d characters.\n", MAXFILELEN-10);
                exit(1);
            }
        } else if (strcmp(argv[arg], "-w") == 0){
            print_eig = 1;
        } else if (strcmp(argv[arg], "-n") == 0){
            if (argv[arg + 1] == NULL || argv[arg + 1][0] == '-') {
                continue;
            }
            arg++;
            errno = 0;
            nperm_long = strtol(argv[arg], &end, 10);
            if(argv[arg] == end || errno != 0) {
                printf("ERROR: Cannot recognize the specified number of replicates to generate!\n");
                exit(1);
            }
            if(nperm_long > NPERM_MAX || nperm_long < 0) {
                printf("ERROR: Number of replicates must be positive and less than %d.\n", NPERM_MAX);
                exit(1);
            } else n_perm = (int) nperm_long;
            nperm_in = 1;
        } else if (strcmp(argv[arg], "--seed") == 0){
            if (argv[arg + 1] == NULL || argv[arg + 1][0] == '-') {
                continue;
            }
            arg++;
            errno= 0;
            seed_rand = strtol(argv[arg], &end, 10);
            if(argv[arg] == end || errno != 0)
                printf("ERROR: Can't recognize the seed speficied! A default seed will be used.\n");
            else seed_in = 1;
        } else if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0){
            printf("Usage: BRASS -p -e -k -o -w -n --seed                       \n"
                    " -p/--ped: name of phenotype file                          \n"
                    " -e/--eigen: name of eigendecomposition file               \n"
                    " -k/--kin: name of kinship coefficients file               \n"
                    " -o/--out: prefix used for the output files                \n"
                    " -w: to write the eigendecomposition results to a file     \n"
                    " -n: specify the number of trait replicates to simulate    \n"
                    " --seed: random seed used for random number generation     \n");
            exit(1);
        } else {
            printf("ERROR: Unrecognized option \"%s\" is ignored!\n", argv[arg]);
            continue;
        }
    }

    // For phenotype/covariate data
    if (!pfile) {
        printf("Default phenotype and covariate data file: %s\n", pheno_filename);
    } else{
        printf("User specified phenotype and covariate data file: %s\n", pheno_filename);
    }

    // For GRM (or its eigen-decomposition)
    if(efile){
        if(kfile){
            printf("WARNING: Only ONE of the options, -e and -r, can be used. Program has ignored -r. Option -e will be taken. \n");
        }
        printf("Program will use eigen-decomposition results provided in the file: %s\n", eig_filename);
        printf("User should make sure this file has the required format.\n\n");
    } else if(kfile){
        printf("Eigen-decomposition will be performed based on the genetic relationship matrix (GRM) provided in the file: %s\n", kin_filename);
        printf("User should make sure this file has the required format.\n\n");
    } else {
        printf("ERROR: If an eigen-decomposition file is not provided by -e, the user needs to provide a genetic relationship matrix (flag -r).\n\n");
        exit(1);
    }

    // Seed for simulations
    if (!seed_in) {
        seed_rand = (long)time(NULL);
        printf("Default seed: %ld\n", seed_rand);
    } else{
        printf("User specified seed: %ld\n", seed_rand);
    }
    // Number of replicates for simulations
    if (!nperm_in) {
        n_perm = 1000;
        printf("Default number of replicates to simulated: %d\n\n", n_perm);
    } else{
        printf("User specified number of replicates to simulate: %ld\n\n", n_perm);
    }

    // For output
    if(!ofile){
        strcpy(param_filename, "BRASS_param.txt");
        strcpy(yrep_filename, "BRASS_perms.txt");
        strcpy(outeig_filename, "BRASS_eig");
    }else{
        strcpy(param_filename, out_prefix);
        strcat(param_filename, "_param.txt");
        strcpy(yrep_filename, out_prefix);
        strcat(yrep_filename, "_perms.txt");
        strcpy(outeig_filename, out_prefix);
        strcat(outeig_filename, "_eig");
    }

    int n_fam = 1, n_person, n_cov, i, j, fam = 1, index_tot, n_total;
    double xi_hat;
    long seed = (-1)*seed_rand;
    srand(seed_rand);

    struct FAMILY *family;
    struct DATA_STRUCT data_struct;
    // Hash table for individual ID
    struct HASH hash;
    hash.ind2ind = kh_init(str);

    // Load phenotype/covariates from file
    data_struct.n_fam = n_fam;
    data_struct.n_perm = n_perm;
    data_struct.tol = NUMTOL;
    printf("Loading phenotype and covariate data...\n");
    format_checks(pheno_filename, &data_struct);
    n_cov =  data_struct.n_cov;
    n_total = data_struct.n_total;

    family = (struct FAMILY *)malloc(sizeof(struct FAMILY) * (n_fam + 1));
    readpheno(pheno_filename, family, &hash, data_struct); // Reads pedigree info
    printf("1 phenotype and %d covariates (intercept not included) found for %d individuals.\n\n", data_struct.n_cov - 1, data_struct.n_total);

    // Get eigendecomposition of GRM
    if(efile){ // If eigendecomposition is provided
        printf("Loading eigen-decomposition results...\n");
        readeigdecomp(eig_filename, family);
        printf("Completed.\n\n");
    } else if (kfile){ // If GRM is provided
        printf("Loading GRM...\n");
        readkin(kin_filename, family, hash);
        printf("GRM loaded.\n");

        printf("Performing eigen-decomposition...\n");
        if(lapa){
            leigcomp_kin(family);
        } else{
            eigcomp_kin(family);
        }
        free_dmatrix((family+fam)->phi, 1, n_total, 1, n_total);

        // Write results to file
        if(print_eig) {
            write_eigendecomp(outeig_filename, family);
            printf("Eigen-decomposition completed. Results stored in: %s\n\n", outeig_filename);
        } else printf("Eigen-decomposition completed.\n\n");
    } else {
        printf("ERROR: Either the GRM or its eigen-decomposition results must be provided!");
        exit(1);
    }

    // Check non-negativeness of eigenvalues
    for (i = 1; i < n_total; i++){
        if((family+fam)->eigvalue[i]< -NUMTOL){
            printf("WARNING: GRM is not symmetric positive semi-definite!\n");
        }
    }

    double time_out;
    clock_t time1, time2, time3, time4;
    //time1 = clock();

    // Fit the null model of CARAT for observed trait
    printf("Estimating parameters under the null...\n");
    xi_hat = compute_gee_res(0, family, data_struct); // Method 0 corresponds to observed trait
    printf("Parameter estimation completed.");

    if(n_perm > 0) {
        // Compute 2nd-order exchangeable vector and transformation matrix to add back the covariance structure to the permuted vector
        pre_permute(1, xi_hat, n_cov, family, param_filename);

        // Starting generation of replicates
        printf("Sampling %d trait replicates...\n", n_perm);
        permute(family, hash, data_struct, yrep_filename);
        printf("Sampling of trait replicates completed. Results stored in %s\n\n", yrep_filename);
    } else {
        pre_permute(0, xi_hat, n_cov, family, param_filename);
        printf("User specified no replicates to be simulated.\n");
    }

    //	/////////////////////////////////////////////////////
    //	////////////////// Free memory //////////////////////
    //	/////////////////////////////////////////////////////
    for (i = 1; i <= n_fam; i++) {
        n_person = family[i].n_total;
        free_dvector((family+i)->trait_obs, 1, n_person);
        free_dvector((family+i)->Z, 1, n_person);
        free_dmatrix((family+i)->cov, 1, n_person, 1, n_cov);
        free_dvector((family+i)->eigvalue, 1, n_person);
        free_dmatrix((family+i)->eigvector, 1, n_person, 1, n_person);
        free_dvector((family+i)->trait, 1, n_person);
        free_dvector((family+i)->Dinv, 1, n_person);
        free_dvector((family+i)->beta_hat, 1, n_cov);
        free_dvector((family+i)->mu_hat_QL, 1, n_person);
        free_dmatrix((family+i)->pre_mult_mat, 1, n_total, 1, n_total - n_cov);
        free_dvector((family+i)->delta, 1, n_total - n_cov);
        if(n_perm > 0) free_dvector((family+i)->trait_rep, 1, n_person);
        for (j = 1; j <= n_person; j++) free((family+i)->indiv[j]);
        free((family+i)->indiv);
    }

    free(family);

    for (int k = kh_begin(hash.ind2ind); k != kh_end(hash.ind2ind); k++){
        if (kh_exist(hash.ind2ind,k)) {
            free((char*) kh_key(hash.ind2ind,k)); /* cast away constness */
        }
    }
    kh_destroy(str,hash.ind2ind);

    return 0;
} /* end main */

void leigcomp_kin(struct FAMILY *family){

    int i, j, fam = 1;
    int n_total = (family + fam)->n_total;

    /* Create a dense matrix for EVD and initalize with GRM entries */
    double *A = dvector(0, n_total * n_total -1);
    for(j = 1 ; j <= n_total; j++)
        for(i = 1 ; i <= n_total ; i++)
            A[(i-1) + (j-1) * n_total] = (family + fam)->phi[i][j];

    (family + fam)->eigvalue = dvector(1, n_total);
    (family + fam)->eigvector = dmatrix(1, n_total, 1, n_total);

    // Initialize arrays for EVD
    int il = 1, evf, info, lwork, liwork, iwkopt, *iwork, *isuppz;
    double wkopt, *work, vl=0, vu=0, *w, *z, abstol=-1.0;
    w = dvector(0, n_total - 1);
    z = dvector(0, n_total * n_total - 1);
    isuppz = ivector(0, 2 * n_total - 1);

    // Get optimal size of the WORK/IWORK arrays
    lwork=-1;
    liwork=-1;
    dsyevr_( "V", "A", "U", &n_total, A, &n_total, &vl, &vu, &il, &n_total, &abstol, &evf, w, z, &n_total, isuppz, &wkopt, &lwork, &iwkopt, &liwork, &info );
    lwork = (int)wkopt;
    liwork = (int)iwkopt;

    // EVD
    work = dvector(0, lwork - 1);
    iwork = ivector(0, liwork - 1);
    dsyevr_( "V", "A", "U", &n_total, A, &n_total, &vl, &vu, &il, &n_total, &abstol, &evf, w, z, &n_total, isuppz, work, &lwork, iwork, &liwork, &info );
    if( info > 0 ) {
        printf( "The algorithm failed to compute eigenvalues for the GRM matrix.\n");
        return;
    }

    // Save the EVD results
    for (j = 1; j <= n_total; j++){
        (family + fam)->eigvalue[j] = w[j - 1];
        for(i = 1 ; i <= n_total ; i++)
            (family + fam)->eigvector[j][i] = z[(j-1)+(i-1)*n_total];
    }

    free_dvector(A, 0, n_total * n_total -1);
    free_dvector(w, 0, n_total -1);
    free_dvector(z, 0, n_total * n_total -1);
    free_ivector(isuppz, 0, 2 * n_total -1);
    free_dvector(work, 0, lwork -1);
    free_ivector(iwork, 0, liwork -1);

}

void eigcomp_kin(struct FAMILY *family){

    int fam = 1;
    int n_total = (family+fam)->n_total;

    // Initialize arrays for EVD result
    (family + fam)->eigvalue = dvector(1, n_total);
    (family + fam)->eigvector = dmatrix(1, n_total, 1, n_total);

    // Map arrays
    MapRowXd Phi_map (&(family+fam)->phi[1][1], n_total, n_total); // GRM
    MapArXd Val_map (&(family + fam)->eigvalue[1], n_total, 1);
    MapRowXd Vec_map (&(family + fam)->eigvector[1][1], n_total, n_total);

    SelfAdjointEigenSolver<MatrixXd> es(Phi_map);
    if (es.info() != Success) {
        printf( "The algorithm failed to compute eigenvalues for the GRM matrix.\n");
        return;
    }

    // Store EVD results
    Val_map = es.eigenvalues().array();
    Vec_map = es.eigenvectors();

}

void write_eigendecomp(char *outeig_filename, struct FAMILY *family){

    int j, k, fam = 1;
    int n_total = (family+fam)->n_total;

    FILE *outeigfile;
    outeigfile = fopen(outeig_filename, "wb");

    for(j = 1; j <= n_total; j++)
        fwrite(&(family + fam)->eigvalue[j], sizeof(double), 1, outeigfile);
    for(j = 1; j <= n_total; j++) {
        for(k = 1; k <= n_total; k++){
            fwrite(&(family + fam)->eigvector[j][k], sizeof(double), 1, outeigfile);
        }
    }
    fclose(outeigfile);
}


// Computes estimated covariate effects assuming Yi's are independent (standard logistic model)
void compute_logistic_coeff(double *yaffec, double **Cov, int totpc, int n_cov, double *beta_init) {
    int i, j, k, itermax = 100, niter = 0;
    double diffpi = 1, tol = NUMTOL;

    MapArXd betaLnew (&beta_init[1], n_cov, 1);
    MapArXd ymap (&yaffec[1], totpc, 1);
    MapRowXd Covmap (&Cov[1][1], totpc, n_cov);

    ArrayXd betaLold(n_cov);
    ArrayXd eta(totpc);
    ArrayXd etaold(totpc);
    ArrayXd piold(totpc);
    ArrayXd pinew(totpc);
    ArrayXd zmap(totpc);
    ArrayXd wmap(totpc);
    ArrayXd XtWZ(n_cov);
    MatrixXd XtW(n_cov, totpc);

    betaLnew.setZero();
    betaLnew(0) = (.5 + ymap.sum()) / (totpc + 1.0);
    betaLnew(0) = log( betaLnew(0) / (1 - betaLnew(0)) );

    eta = (Covmap * betaLnew.matrix()).array();
    pinew = (1 - 1/(eta.exp() + 1));

    while (niter++ < itermax) {
        wmap = (pinew*(1-pinew));
        zmap = eta + (ymap - pinew)/wmap;
        XtW = Covmap.transpose() * wmap.matrix().asDiagonal();
        betaLold = ((XtW * Covmap).colPivHouseholderQr().solve(XtW * zmap.matrix()).array()).array();
        etaold = (Covmap * betaLold.matrix()).array();
        if( (betaLnew-betaLold).square().sum() < tol) break;
        piold = (1 - 1/(etaold.exp()+1));
        if( (pinew - piold).square().sum() < tol) break;

        betaLnew = betaLold;
        eta = etaold;
        pinew = piold;
    }
    //cout << betaLnew << endl << endl;
    // Done with logistic
}

void estimate_gee_parameters(int method, int totpc, double **V, double *L, double *Dinv, struct FAMILY *family, int n_cov, int n_fam, double *xi, double *betahat) {

    int i, j, k, fam = 1;
    int n_total = totpc;
    double *M, f_lo;

    param_func_xi parameters = { n_total, V, L, Dinv, family, betahat, n_fam, n_cov };

    function_xi *func_xi = (function_xi *)malloc(sizeof(function_xi));
    func_xi->parameters = parameters;
    func_xi->func_of_xi = &eqn3_lhs_minus_rhs;

    brent(func_xi, xi);

    MapArXd Lmap(&L[1], n_total, 1);
    MapArXd Dinv_map(&Dinv[1], n_total, 1);
    Dinv_map = (*xi * Lmap + 1 - *xi).inverse();
    compute_gee_betahat(V, Dinv, betahat, family, n_fam, n_cov);

    if(method == 0) { 
        // Save mu for observed trait
        MapRowXd x_map(&(family + fam)->cov[1][1], n_total, n_cov);
        MapRowXd beta(&betahat[1], n_cov, 1);
        MapArXd mu(&(family + fam)->mu_hat_QL[1], n_total, 1);
        mu = 1 - 1 / ((x_map * beta).array().exp() + 1);
    }

    free(func_xi);
}

// totpc is n_total, xi is heritability
double target_func_xi(double xi, int totpc, double **V, double *L, double *Dinv, struct FAMILY *family, double *current_beta, int n_fam, int n_cov) {

    int j, fam, flag_betahat = -1;
    double fvalue;

    for(fam = 1; fam <= n_fam; fam++) {
        MapArXd Lmap(&L[1], (family + fam)->n_total, 1);
        MapArXd Dinv_map(&Dinv[1], (family + fam)->n_total, 1);

        Dinv_map = (xi * Lmap + 1 - xi).inverse();
    }

    flag_betahat = compute_gee_betahat(V, Dinv, current_beta, family, n_fam, n_cov);

    if (flag_betahat == 0) {
        fvalue = MISSVAL_FUNC_XI;
    }
    else
        // Refers back to eqn 10 on sheng paper
        //compute rhs_eq3 = trace(Sigma_inverse * (Phi - I)) = tr(Dinv*(L-1))
        //compute_gee_betahat => Z = Dinv * Vt *((Y-u)/M)
        //compute lhs_eq3 = Zt * (L-1) * Z = sum( diag(L-I)*Z^2)
        // fvalue = lhs - rhs
        for (fvalue = 0, fam = 1; fam <= n_fam; fam++) {
            MapArXd Zmap ( &((family+fam)->Z)[1], (family+fam)->n_total, 1);
            MapArXd Dinv_map ( &(Dinv[1]), (family+fam)->n_total, 1);
            MapArXd Lmap ( &L[1], (family+fam)->n_total, 1);

            fvalue += (Zmap.array().square() * (Lmap - 1).array()).sum() - (Dinv_map.array() * (Lmap - 1).array()).sum();
        }

    return fabs(fvalue);
}

int compute_gee_betahat(double **V, double *Dinv, double *current_beta, struct FAMILY *family, int n_fam, int n_cov) {

    int i, j, k, l, fam;
    double betadiff=1, tmp, tol = 1e-8;
    // printf("Current xi is: %f.\n",xi);
    int niter = 0, flag_while = 1;
    int n_person = (family+1)->n_total;

    MatrixXd beta (n_cov, 1);
    MatrixXd XtWX (n_cov, n_cov);
    MatrixXd XtM_SigmaInv_MInv_y(n_cov, 1);
    MatrixXd beta_new (n_cov, 1);

    while (betadiff > tol) { // Use difference in beta as criterion
        niter++;
        if (niter > 50) {
            flag_while = 0;
            break;
        }
        for (i = 1; i <= n_cov; i++) beta(i - 1, 0) = current_beta[i];

        XtWX.setZero();
        XtM_SigmaInv_MInv_y.setZero();

        //compute mu_f
        for (fam = 1; fam <= n_fam; fam++) {
            MapArXd y_map ( &((family+fam)->trait)[1], n_person, 1);
            MapRowXd x_map (&((family+fam)->cov)[1][1], n_person, n_cov);
            MapArXd Zmap ( &((family+fam)->Z)[1], n_person, 1);

            MapRowXd D_inv_map (&(Dinv[1]), n_person, 1);
            MapColXd V_map (&V[1][1], n_person, n_person); // Vt

            ArrayXd mu = (1 - 1/((x_map * beta).array().exp() + 1));
            ArrayXd M = (mu * (1-mu)).sqrt();

            MatrixXd UtMx = V_map * M.matrix().asDiagonal() * x_map;
            XtWX += UtMx.transpose() * D_inv_map.asDiagonal() * UtMx;
            Zmap = D_inv_map.array() * (V_map * ((y_map - mu) / M).matrix()).array();
            XtM_SigmaInv_MInv_y += UtMx.transpose() * Zmap.matrix();
        }

        //solve quasi-likelihood equations for beta using Newton's method with Fisher scoring
        beta_new = beta + XtWX.colPivHouseholderQr().solve(XtM_SigmaInv_MInv_y);
        betadiff = (beta - beta_new).array().abs().sum();

        for (i = 0; i <n_cov; i++) current_beta[i+1] = beta_new(i,0);
    }

    for (i = 0; i < n_cov; i++) {
        if (std::isnan(beta_new(i, 0))){
            flag_while = 0;
            current_beta[i + 1] = 0;
        }
    }
    return flag_while;
}

// For permutation-based replicates of CARAT
void pre_permute(int get_perm, double xi, int n_cov, struct FAMILY *family, char *filename) {

    int i, fam = 1;
    int n_total = (family+fam)->n_total;
    double *betahat = (family + fam)->beta_hat;

    MapArXd Ymap (&(family + fam)->trait[1], n_total, 1);
    MapRowXd x_map (&(family+fam)->cov[1][1], n_total, n_cov);

    MapArXd eval (&(family+fam)->eigvalue[1], n_total, 1);
    MapRowXd evec (&(family+fam)->eigvector[1][1], n_total, n_total);

    MapArXd mu(&(family + fam)->mu_hat_QL[1], n_total, 1);
    MapRowXd mapV1 (&(family+fam)->pre_mult_mat[1][1], n_total, n_total - n_cov);
    MapArXd delta(&(family + fam)->delta[1], n_total - n_cov, 1);

    // Concatenate (M*X, (Y-mu))
    ArrayXd M = (mu * (1 - mu)).sqrt();
    MatrixXd MX_Yres (x_map.rows(), x_map.cols() + 1);
    MX_Yres << M.matrix().asDiagonal() * x_map, ((Ymap - mu)/M).matrix() ;

    // Pre-multiply by C^{-t} where Sigma = (xi * phi + 1 - xi) = ULUt= C^t * C
    ArrayXd Lsqrt = (xi * eval + 1 - xi).sqrt();
    MatrixXd cholSigInv_Mx_Yres = evec * (1/Lsqrt).matrix().asDiagonal() * evec.transpose() * MX_Yres;

    // Take svd of W = C^{-t} * M * X
    JacobiSVD<MatrixXd> svdW( cholSigInv_Mx_Yres.leftCols(x_map.cols()), ComputeFullU | ComputeFullV);

    // Print estimates + SE for parameters
    ArrayXd SEs =  (svdW.matrixV() * (1/svdW.singularValues().array().square()).matrix().asDiagonal() * svdW.matrixV().transpose()).diagonal().array().sqrt();
    FILE *ests;
    ests = fopen(filename, "w");

    fprintf(ests, "Parameter\tNull_Estimate\tSE\n");
    fprintf(ests, "Variance_Parameter_xi\t%.17g\tNA\n", xi);
    fprintf(ests, "Intercept\t%.17g\t%.17g\n", betahat[1], SEs(0));
    for (i = 2; i <= n_cov; i++) fprintf(ests, "Covariate_%d\t%.17g\t%.17g\n", i - 1, betahat[i], SEs(i-1));

    fclose(ests);
    printf(" Results stored in %s\n\n", filename);
    ////////////////////////////////////////

    // Arrays used to generate trait replicates
    if(get_perm == 1) {
        printf("Computing uncorrelated residuals and un-scaling matrix...\n");
        //// Matrices and vectors needed for un-centering and un-scaling
        mapV1 = svdW.matrixU().rightCols(n_total - n_cov);
        delta = (mapV1.transpose() * cholSigInv_Mx_Yres.rightCols(1)).array();
        mapV1 = M.matrix().asDiagonal() * evec * Lsqrt.matrix().asDiagonal() * evec.transpose() * mapV1;
    }

}


void get_quant_rep(struct FAMILY *family, struct DATA_STRUCT data_struct){

    int fam = 1;
    int n_total = data_struct.n_total;
    int n_cov = data_struct.n_cov;
    int nonzero = n_total - n_cov;
    double *delta = (family + fam)->delta;

    // Permute entries of delta
    naive_shuffle(delta, nonzero);

    // Arrays needed
    MapArXd Ymap (&(family+fam)->trait_rep[1], n_total, 1);
    MapArXd mu(&(family + fam)->mu_hat_QL[1], n_total, 1);
    MapRowXd mapV1 (&(family+fam)->pre_mult_mat[1][1], n_total, nonzero);
    MapArXd mapDelta(&(family + fam)->delta[1], nonzero, 1);

    Ymap = mu +  (mapV1 * mapDelta.matrix()).array();

}


// Assume starting index is 1
void naive_shuffle(double *trait, int n_total){
    random_shuffle(&trait[1], &trait[1 + n_total]);
}
